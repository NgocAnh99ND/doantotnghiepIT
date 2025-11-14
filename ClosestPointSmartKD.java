import java.util.*;
import java.util.function.Supplier;

public class ClosestPointSmartKD {

    // ========================= POINT =========================
    static class Point {
        double lat, lon;
        Point(double lat, double lon) { this.lat = lat; this.lon = lon; }
    }

    // ========================= RECT =========================
    static class Rect {
        double minLat, minLon, maxLat, maxLon;

        Rect(double minLat, double minLon, double maxLat, double maxLon) {
            this.minLat = minLat;
            this.minLon = minLon;
            this.maxLat = maxLat;
            this.maxLon = maxLon;
        }

        double width()  { return maxLon - minLon; }
        double height() { return maxLat - minLat; }

        Point center() {
            return new Point((minLat + maxLat)/2, (minLon + maxLon)/2);
        }

        double diagonalMeters() {
            return haversineMeters(
                    new Point(minLat, minLon),
                    new Point(maxLat, maxLon)
            );
        }

        boolean intersects(Point A, Point B) {
            return !(Math.max(A.lat, B.lat) < minLat ||
                     Math.min(A.lat, B.lat) > maxLat ||
                     Math.max(A.lon, B.lon) < minLon ||
                     Math.min(A.lon, B.lon) > maxLon);
        }
    }

    // ========================= KD-TREE =========================
    static class KDNode {
        Point p;
        KDNode left, right;
        boolean splitLat;
        KDNode(Point p, boolean splitLat) { this.p = p; this.splitLat = splitLat; }
    }

    static KDNode buildKDTree(List<Point> pts, boolean splitLat) {
        if (pts.isEmpty()) return null;

        pts.sort((a, b) -> splitLat
                ? Double.compare(a.lat, b.lat)
                : Double.compare(a.lon, b.lon));

        int mid = pts.size()/2;
        KDNode node = new KDNode(pts.get(mid), splitLat);

        node.left = buildKDTree(pts.subList(0, mid), !splitLat);
        node.right = buildKDTree(pts.subList(mid+1, pts.size()), !splitLat);

        return node;
    }

    static Point kdNearest(KDNode node, Point target, Point best, double[] bestDist) {
        if (node == null) return best;

        double d = haversineMeters(target, node.p);
        if (d < bestDist[0]) {
            bestDist[0] = d;
            best = node.p;
        }

        boolean goLeft = node.splitLat
                ? target.lat < node.p.lat
                : target.lon < node.p.lon;

        KDNode first  = goLeft ? node.left  : node.right;
        KDNode second = goLeft ? node.right : node.left;

        best = kdNearest(first, target, best, bestDist);

        double axisDist = node.splitLat
                ? Math.abs(target.lat - node.p.lat)*111000
                : Math.abs(target.lon - node.p.lon)*111000;

        if (axisDist < bestDist[0])
            best = kdNearest(second, target, best, bestDist);

        return best;
    }

    // ========================= HAVERSINE =========================
    static double haversineMeters(Point a, Point b) {
        double R = 6371000.0;
        double dLat = Math.toRadians(b.lat - a.lat);
        double dLon = Math.toRadians(b.lon - a.lon);
        double lat1 = Math.toRadians(a.lat);
        double lat2 = Math.toRadians(b.lat);

        double h = Math.sin(dLat/2)*Math.sin(dLat/2)
                + Math.cos(lat1)*Math.cos(lat2)
                * Math.sin(dLon/2)*Math.sin(dLon/2);

        return 2 * R * Math.asin(Math.sqrt(h));
    }

    // ========================= NỘI SUY =========================
    static Point interpolate(Point A, Point B, double f) {
        return new Point(A.lat + f*(B.lat - A.lat),
                         A.lon + f*(B.lon - A.lon));
    }

    static List<Point> adaptiveResample(List<Point> src, double spacingM) {
        List<Point> out = new ArrayList<>();
        out.add(src.get(0));
        for (int i=0; i<src.size()-1; i++) {
            Point A = src.get(i);
            Point B = src.get(i+1);
            double d = haversineMeters(A, B);

            if (d > spacingM) {
                int n = (int)(d/spacingM);
                for (int j=1;j<=n;j++)
                    out.add(interpolate(A, B, (double)j/(n+1)));
            }
            out.add(B);
        }
        return out;
    }

    // ========================= TẠO TUYẾN 200 ĐIỂM – CONG – 100KM =========================
    static List<Point> generateRoute200() {
        List<Point> pts = new ArrayList<>();

        double startLat = 21.028;
        double startLon = 105.778;

        double totalKm = 100.0;
        double latKmToDeg = 1.0 / 111.0;
        double lonKmToDeg = 1.0 / (111.0 * Math.cos(Math.toRadians(startLat)));

        for (int i = 0; i < 200; i++) {
            double t = (double)i / 199.0;

            double dist = t * totalKm;

            double lat = startLat - dist * latKmToDeg
                    + 0.01 * Math.sin(t * Math.PI * 6)
                    + 0.005 * Math.sin(t * Math.PI * 12);

            double lon = startLon + dist * lonKmToDeg
                    + 0.015 * Math.cos(t * Math.PI * 4)
                    + 0.005 * Math.sin(t * Math.PI * 8);

            pts.add(new Point(lat, lon));
        }

        return pts;
    }

    // ========================= SHRINK RECT =========================
    static Rect shrinkRect(Rect rect, Point customer, double stopMeters) {
        Rect cur = rect;
        while (cur.diagonalMeters() > stopMeters) {

            boolean splitHoriz = (cur.width() >= cur.height());

            Rect left, right;

            if (splitHoriz) {
                double mid = (cur.minLon + cur.maxLon)/2;
                left  = new Rect(cur.minLat, cur.minLon, cur.maxLat, mid);
                right = new Rect(cur.minLat, mid, cur.maxLat, cur.maxLon);
            } else {
                double mid = (cur.minLat + cur.maxLat)/2;
                left  = new Rect(cur.minLat, cur.minLon, mid, cur.maxLon);
                right = new Rect(mid, cur.minLon, cur.maxLat, cur.maxLon);
            }

            Point cL = left.center();
            Point cR = right.center();

            double dL = haversineMeters(customer, cL);
            double dR = haversineMeters(customer, cR);

            Rect chosen = (dL < dR) ? left : right;
            Point centerChosen = (dL < dR) ? cL : cR;

            double axisDist = splitHoriz ?
                    Math.abs(centerChosen.lon - (cur.minLon + cur.maxLon)/2)*111000
                    :
                    Math.abs(centerChosen.lat - (cur.minLat + cur.maxLat)/2)*111000;

            double ratio = Math.min(dL, dR) / axisDist;

            cur = chosen;
        }

        return cur;
    }

    // ========================= FIND CLOSEST =========================
    static Point findSmartClosest(Point customer, List<Point> route, double stopMeters, double spacingM) {

        double minLat=1e9, minLon=1e9, maxLat=-1e9, maxLon=-1e9;
        for (Point p: route){
            minLat=Math.min(minLat,p.lat);
            minLon=Math.min(minLon,p.lon);
            maxLat=Math.max(maxLat,p.lat);
            maxLon=Math.max(maxLon,p.lon);
        }

        Rect rect = new Rect(minLat,minLon,maxLat,maxLon);

        Rect small = shrinkRect(rect, customer, stopMeters);

        List<Point> pruned = new ArrayList<>();
        for (int i=0;i<route.size()-1;i++){
            if (small.intersects(route.get(i),route.get(i+1))) {
                pruned.add(route.get(i));
                pruned.add(route.get(i+1));
            }
        }
        if (pruned.isEmpty()) pruned = route;

        List<Point> dense = adaptiveResample(pruned, spacingM);

        KDNode root = buildKDTree(dense,true);
        double[] bestDist = {Double.MAX_VALUE};
        return kdNearest(root,customer,null,bestDist);
    }

    // ========================= TIMER =========================
    static <T> T measureTime(String label, Supplier<T> task){
        long s=System.nanoTime();
        T r=task.get();
        long e=System.nanoTime();
        System.out.printf("⏰ %s: %.3f ms%n",label,(e-s)/1e6);
        return r;
    }

    // ========================= MAIN =========================
    public static void main(String[] args) {

        List<Point> route = generateRoute200();
        System.out.println("📌 Tổng số điểm lộ trình: " + route.size());

        Point customer = new Point(20.992, 105.816);

        Point result = measureTime("Smart KD-tree nearest", () ->
            findSmartClosest(customer, route, 500.0, 50.0)
        );

        double dist = haversineMeters(customer, result);

        System.out.println("\n===== KẾT QUẢ =====");
        System.out.printf("📍 Customer: (%.6f, %.6f)%n", customer.lat, customer.lon);
        System.out.printf("🎯 Closest:  (%.6f, %.6f)%n", result.lat, result.lon);
        System.out.printf("📏 Distance: %.2f m%n", dist);
    }
}
