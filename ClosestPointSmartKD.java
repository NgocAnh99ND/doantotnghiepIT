import java.util.*;

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
        Point center()  { return new Point((minLat+maxLat)/2, (minLon+maxLon)/2); }

        double diagonalMeters() {
            return ClosestPointSmartKD.haversineMeters(
                    new Point(minLat,minLon),
                    new Point(maxLat,maxLon)
            );
        }

        boolean intersects(Point A, Point B) {
            return !(Math.max(A.lat,B.lat) < minLat ||
                    Math.min(A.lat,B.lat) > maxLat ||
                    Math.max(A.lon,B.lon) < minLon ||
                    Math.min(A.lon,B.lon) > maxLon);
        }
    }

    // ========================= KD-TREE =========================
    static class KDNode {
        Point p;
        KDNode left, right;
        boolean splitLat;
        KDNode(Point p, boolean splitLat){ this.p=p; this.splitLat=splitLat; }
    }

    static KDNode buildKDTree(List<Point> pts, boolean splitLat){
        if(pts.isEmpty()) return null;

        pts.sort((a,b)-> splitLat
                ? Double.compare(a.lat,b.lat)
                : Double.compare(a.lon,b.lon));

        int mid = pts.size()/2;
        KDNode node = new KDNode(pts.get(mid), splitLat);

        node.left  = buildKDTree(pts.subList(0,mid), !splitLat);
        node.right = buildKDTree(pts.subList(mid+1,pts.size()), !splitLat);

        return node;
    }

    static Point kdNearest(KDNode node, Point target, Point best, double[] bestDist){
        if(node == null) return best;

        double d = haversineMeters(target, node.p);
        if(d < bestDist[0]){
            bestDist[0] = d;
            best = node.p;
        }

        boolean goLeft = node.splitLat
                ? target.lat < node.p.lat
                : target.lon < node.p.lon;

        KDNode first  = goLeft ? node.left  : node.right;
        KDNode second = goLeft ? node.right : node.left;

        best = kdNearest(first, target, best, bestDist);

        double axisDist = node.splitLat ?
                Math.abs(target.lat - node.p.lat)*111000 :
                Math.abs(target.lon - node.p.lon)*111000;

        if(axisDist < bestDist[0])
            best = kdNearest(second, target, best, bestDist);

        return best;
    }

    // ========================= HAVERSINE =========================
    static double haversineMeters(Point a, Point b){
        double R=6371000;
        double dLat=Math.toRadians(b.lat - a.lat);
        double dLon=Math.toRadians(b.lon - a.lon);
        double lat1=Math.toRadians(a.lat);
        double lat2=Math.toRadians(b.lat);

        double h = Math.sin(dLat/2)*Math.sin(dLat/2)
                + Math.cos(lat1)*Math.cos(lat2)
                * Math.sin(dLon/2)*Math.sin(dLon/2);

        return 2*R*Math.asin(Math.sqrt(h));
    }

    // ========================= NỘI SUY =========================
    static Point interpolate(Point A, Point B, double f){
        return new Point(
                A.lat + f*(B.lat-A.lat),
                A.lon + f*(B.lon-A.lon)
        );
    }

    static List<Point> adaptiveResample(List<Point> src, double spacingM){
        List<Point> out=new ArrayList<>();
        if (src.isEmpty()) return out;

        out.add(src.get(0));

        for(int i=0;i<src.size()-1;i++){
            Point A=src.get(i), B=src.get(i+1);
            double d = haversineMeters(A,B);

            if(d > spacingM){
                int n = (int)(d/spacingM);
                for(int j=1;j<=n;j++)
                    out.add(interpolate(A,B,(double)j/(n+1)));
            }
            out.add(B);
        }
        return out;
    }

    // ========================= TẠO TUYẾN 10 ĐIỂM =========================
    static List<Point> generateRoute10(){
        return Arrays.asList(
            new Point(21.028, 105.778),
            new Point(21.025, 105.785),
            new Point(21.021, 105.795),
            new Point(21.018, 105.805),
            new Point(21.015, 105.815),
            new Point(21.011, 105.825),
            new Point(21.007, 105.835),
            new Point(21.004, 105.845),
            new Point(21.000, 105.855),
            new Point(20.994, 105.866)
        );
    }

    // ========================= MAIN =========================
    public static void main(String[] args) {

        // 1) Lộ trình 10 điểm
        List<Point> route = generateRoute10();

        // 2) Khách
        Point customer = new Point(20.992,105.816);

        // 3) Nội suy trên toàn tuyến
        List<Point> dense = adaptiveResample(route, 50.0);

        // 4) BUILD KD-tree
        long t1 = System.nanoTime();
        KDNode root = buildKDTree(dense, true);
        long t2 = System.nanoTime();
        System.out.printf("Thoi gian build KD-tree: %.3f ms%n", (t2-t1)/1e6);

        // 5) TÌM GẦN NHẤT
        long s = System.nanoTime();
        double[] bestDist = {Double.MAX_VALUE};
        Point closest = kdNearest(root, customer, null, bestDist);
        long e = System.nanoTime();
        System.out.printf("Thoi gian chay kdNearest: %.3f ms%n", (e-s)/1e6);

        System.out.println("\n===== KET QUA =====");
        System.out.printf("Customer: (%.6f, %.6f)%n", customer.lat, customer.lon);
        System.out.printf("Closest : (%.6f, %.6f)%n", closest.lat, closest.lon);
        System.out.printf("Distance: %.2f m%n", bestDist[0]);
    }
}
