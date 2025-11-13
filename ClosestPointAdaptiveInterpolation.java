import java.util.*;
import java.util.function.Supplier;

public class ClosestPointAdaptiveInterpolation {

    static class Point {
        double lat, lon;
        Point(double lat, double lon) { this.lat = lat; this.lon = lon; }
    }

    static class Rect {
        double minLat, minLon, maxLat, maxLon;
        Rect(double minLat, double minLon, double maxLat, double maxLon) {
            this.minLat = minLat; this.minLon = minLon;
            this.maxLat = maxLat; this.maxLon = maxLon;
        }
        double width()  { return maxLon - minLon; }
        double height() { return maxLat - minLat; }
        Point center()  { return new Point((minLat + maxLat)/2, (minLon + maxLon)/2); }
        double diagonalMeters() {
            return haversineMeters(new Point(minLat, minLon), new Point(maxLat, maxLon));
        }
    }

    // ===== HAVERSINE DISTANCE (m) =====
    static double haversineMeters(Point a, Point b) {
        double R = 6371000.0; // Bán kính Trái đất (m)
        double dLat = Math.toRadians(b.lat - a.lat);
        double dLon = Math.toRadians(b.lon - a.lon);
        double lat1 = Math.toRadians(a.lat);
        double lat2 = Math.toRadians(b.lat);
        double h = Math.pow(Math.sin(dLat / 2), 2)
                + Math.pow(Math.sin(dLon / 2), 2) * Math.cos(lat1) * Math.cos(lat2);
        return 2 * R * Math.asin(Math.sqrt(h));
    }

    // ===== CHIA VÙNG THEO CẠNH DÀI HƠN =====
    static List<Rect> splitByLongestSide(Rect r) {
        if (r.width() >= r.height()) {
            double midLon = (r.minLon + r.maxLon) / 2;
            return Arrays.asList(
                    new Rect(r.minLat, r.minLon, r.maxLat, midLon),
                    new Rect(r.minLat, midLon, r.maxLat, r.maxLon)
            );
        } else {
            double midLat = (r.minLat + r.maxLat) / 2;
            return Arrays.asList(
                    new Rect(r.minLat, r.minLon, midLat, r.maxLon),
                    new Rect(midLat, r.minLon, r.maxLat, r.maxLon)
            );
        }
    }

    // ===== TÍNH ĐIỂM NỘI SUY =====
    static Point interpolate(Point A, Point B, double f) {
        return new Point(A.lat + f * (B.lat - A.lat),
                         A.lon + f * (B.lon - A.lon));
    }

    // ===== CHỈ NỘI SUY NẾU ĐOẠN > spacingM =====
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

    // ===== TÌM ĐIỂM GẦN P NHẤT TRÊN ĐOẠN AB =====
    static Point closestPointOnSegment(Point A, Point B, Point P) {
        double x1 = A.lon, y1 = A.lat, x2 = B.lon, y2 = B.lat;
        double x0 = P.lon, y0 = P.lat;
        double dx = x2 - x1, dy = y2 - y1;
        if (dx == 0 && dy == 0) return A;
        double t = ((x0 - x1)*dx + (y0 - y1)*dy) / (dx*dx + dy*dy);
        t = Math.max(0, Math.min(1, t));
        return new Point(y1 + t*dy, x1 + t*dx);
    }

    // ===== CHÍNH: CHIA VÙNG + DUYỆT =====
    static Point findClosestOnRoute(Point customer, List<Point> route, Rect initialRect,
                                    double stopMeters, double spacingM) {
        Rect cur = initialRect;
        int iter = 0;

        // B1 + B2: chia đôi cho đến khi vùng < stopMeters
        while (cur.diagonalMeters() > stopMeters) {
            List<Rect> halves = splitByLongestSide(cur);
            Rect left = halves.get(0);
            Rect right = halves.get(1);

            double dLeft = haversineMeters(customer, left.center());
            double dRight = haversineMeters(customer, right.center());
            cur = (dLeft < dRight) ? left : right;
            iter++;
        }
        System.out.printf("🔁 Số lần chia vùng: %d%n", iter);

        // B3: chỉ nội suy nếu cần
        List<Point> denseRoute = adaptiveResample(route, spacingM);

        // B4: duyệt tuyến để tìm điểm gần khách nhất
        double bestDist = Double.MAX_VALUE;
        Point best = null;
        for (int i = 0; i < denseRoute.size() - 1; i++) {
            Point A = denseRoute.get(i);
            Point B = denseRoute.get(i + 1);
            Point C = closestPointOnSegment(A, B, customer);
            double d = haversineMeters(customer, C);
            if (d < bestDist) {
                bestDist = d;
                best = C;
            }
        }
        System.out.printf("📏 Khoảng cách nhỏ nhất: %.1f m%n", bestDist);
        return best;
    }

    // ===== HÀM ĐO THỜI GIAN =====
    static <T> T measureTime(String label, Supplier<T> task) {
        long start = System.nanoTime();
        T result = task.get();
        long end = System.nanoTime();
        System.out.printf("⏰ %s: %.3f ms%n", label, (end - start)/1_000_000.0);
        return result;
    }

    // ===== MAIN =====
    public static void main(String[] args) {
        // Lộ trình tài xế (A → J)
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

        // Vị trí khách
        Point customer = new Point(20.992, 105.816);

        // Vùng bao quanh tuyến
        Rect routeRect = new Rect(20.994, 105.778, 21.028, 105.866);

        // Thực thi
        Point result = measureTime("Tìm điểm gần khách nhất trên tuyến", () ->
            findClosestOnRoute(customer, route, routeRect, 500.0, 50.0)
        );

        double distM = haversineMeters(customer, result);
        System.out.println("\n=== KẾT QUẢ CUỐI CÙNG ===");
        System.out.printf("📍 Khách     : (%.6f, %.6f)%n", customer.lat, customer.lon);
        System.out.printf("🎯 Gần nhất  : (%.6f, %.6f)%n", result.lat, result.lon);
        System.out.printf("📏 Khoảng cách: %.1f m%n", distM);
        System.out.printf("🔗 Mở trên Google Maps: https://www.google.com/maps?q=%.6f,%.6f%n",
                result.lat, result.lon);
    }
}
