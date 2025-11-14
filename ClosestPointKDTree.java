import java.util.*;
import java.util.function.Supplier;

public class ClosestPointKDTree {

    // ========================= POINT =========================
    static class Point {
        double lat, lon;
        Point(double lat, double lon) { this.lat = lat; this.lon = lon; }
    }

    // ========================= KD-TREE NODE =========================
    static class KDNode {
        Point p;
        KDNode left, right;
        boolean splitLat; // true = split by lat, false = split by lon

        KDNode(Point p, boolean splitLat) {
            this.p = p;
            this.splitLat = splitLat;
        }
    }

    // ========================= KD-TREE BUILD =========================
    static KDNode buildKDTree(List<Point> pts, boolean splitLat) {
        if (pts.isEmpty()) return null;

        pts.sort((a, b) -> splitLat
                ? Double.compare(a.lat, b.lat)
                : Double.compare(a.lon, b.lon));

        int mid = pts.size() / 2;
        KDNode node = new KDNode(pts.get(mid), splitLat);

        node.left  = buildKDTree(pts.subList(0, mid), !splitLat);
        node.right = buildKDTree(pts.subList(mid + 1, pts.size()), !splitLat);

        return node;
    }

    // ========================= KD-TREE QUERY =========================
    static Point kdNearest(KDNode node, Point target, Point best, double[] bestDist) {
        if (node == null) return best;

        double d = haversineMeters(target, node.p);
        if (d < bestDist[0]) {
            bestDist[0] = d;
            best = node.p;
        }

        boolean goLeftFirst;
        if (node.splitLat)
            goLeftFirst = target.lat < node.p.lat;
        else
            goLeftFirst = target.lon < node.p.lon;

        KDNode first  = goLeftFirst ? node.left : node.right;
        KDNode second = goLeftFirst ? node.right : node.left;

        best = kdNearest(first, target, best, bestDist);

        // kiểm tra có cần duyệt nhánh còn lại không
        double axisDist = node.splitLat ?
                Math.abs(target.lat - node.p.lat) * 111000 :
                Math.abs(target.lon - node.p.lon) * 111000;

        if (axisDist < bestDist[0]) {
            best = kdNearest(second, target, best, bestDist);
        }

        return best;
    }

    // ========================= HAVERSINE (m) =========================
    static double haversineMeters(Point a, Point b) {
        double R = 6371000.0;
        double dLat = Math.toRadians(b.lat - a.lat);
        double dLon = Math.toRadians(b.lon - a.lon);
        double lat1 = Math.toRadians(a.lat);
        double lat2 = Math.toRadians(b.lat);

        double h = Math.sin(dLat/2)*Math.sin(dLat/2)
                + Math.cos(lat1)*Math.cos(lat2)*Math.sin(dLon/2)*Math.sin(dLon/2);

        return 2 * R * Math.asin(Math.sqrt(h));
    }

    // ========================= NỘI SUY 50m =========================
    static Point interpolate(Point A, Point B, double f) {
        return new Point(
                A.lat + f * (B.lat - A.lat),
                A.lon + f * (B.lon - A.lon)
        );
    }

    static List<Point> adaptiveResample(List<Point> src, double spacingM) {
        List<Point> out = new ArrayList<>();
        if (src.isEmpty()) return out;
        out.add(src.get(0));

        for (int i = 0; i < src.size() - 1; i++) {
            Point A = src.get(i);
            Point B = src.get(i + 1);
            double d = haversineMeters(A, B);

            if (d > spacingM) {
                int n = (int) Math.floor(d / spacingM);
                for (int j = 1; j <= n; j++) {
                    double f = (double) j / (n + 1);
                    out.add(interpolate(A, B, f));
                }
            }
            out.add(B);
        }
        return out;
    }

    // ========================= MAIN FIND FUNCTION =========================
    static Point findClosestKD(Point customer, List<Point> route, double spacingM) {

        // bước 1: nội suy
        List<Point> dense = adaptiveResample(route, spacingM);

        // bước 2: xây KD-tree
        KDNode root = buildKDTree(dense, true);

        // bước 3: truy vấn nearest
        double[] bestDist = { Double.MAX_VALUE };
        return kdNearest(root, customer, null, bestDist);
    }

    // ========================= TIỆN ÍCH =========================
    static <T> T measureTime(String label, Supplier<T> task) {
        long start = System.nanoTime();
        T result = task.get();
        long end = System.nanoTime();
        System.out.printf("⏰ %s: %.3f ms%n", label, (end - start)/1_000_000.0);
        return result;
    }

    // ========================= DEMO MAIN =========================
    public static void main(String[] args) {

        List<Point> route = Arrays.asList(
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

        Point customer = new Point(20.992, 105.816);

        Point result = measureTime("KD-Tree nearest neighbor", () ->
                findClosestKD(customer, route, 50.0)
        );

        double distM = haversineMeters(customer, result);

        System.out.println("\n===== KẾT QUẢ CUỐI CÙNG =====");
        System.out.printf("📍 Khách     : (%.6f, %.6f)%n", customer.lat, customer.lon);
        System.out.printf("🎯 Gần nhất  : (%.6f, %.6f)%n", result.lat, result.lon);
        System.out.printf("📏 Khoảng cách: %.1f m%n", distM);
    }
}
