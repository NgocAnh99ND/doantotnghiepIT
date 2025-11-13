import java.util.*;
import java.util.function.Supplier;

public class NearbyDriversFinder {

    static class Point {
        double lat, lon;

        Point(double lat, double lon) {
            this.lat = lat;
            this.lon = lon;
        }
    }

    static class Rect {
        double minLat, minLon, maxLat, maxLon;

        Rect(double minLat, double minLon, double maxLat, double maxLon) {
            this.minLat = minLat;
            this.minLon = minLon;
            this.maxLat = maxLat;
            this.maxLon = maxLon;
        }
    }

    // ✅ Hàm tính khoảng cách Haversine (km)
    static double haversine(Point a, Point b) {
        double R = 6371.0; // bán kính Trái Đất (km)
        double dLat = Math.toRadians(b.lat - a.lat);
        double dLon = Math.toRadians(b.lon - a.lon);
        double lat1 = Math.toRadians(a.lat);
        double lat2 = Math.toRadians(b.lat);
        double h = Math.pow(Math.sin(dLat / 2), 2)
                + Math.pow(Math.sin(dLon / 2), 2) * Math.cos(lat1) * Math.cos(lat2);
        return 2 * R * Math.asin(Math.sqrt(h));
    }

    // ✅ Tạo vùng tìm kiếm 1 km quanh khách
    static Rect makeSearchRect(Point center, double radiusKm) {
        double deltaLat = radiusKm / 111.0;
        double deltaLon = radiusKm / (111.0 * Math.cos(Math.toRadians(center.lat)));
        return new Rect(center.lat - deltaLat / 2, center.lon - deltaLon / 2,
                center.lat + deltaLat / 2, center.lon + deltaLon / 2);
    }

    // ✅ Tính bounding box của lộ trình
    static Rect boundingRect(List<Point> route) {
        double minLat = Double.MAX_VALUE, minLon = Double.MAX_VALUE;
        double maxLat = -Double.MAX_VALUE, maxLon = -Double.MAX_VALUE;
        for (Point p : route) {
            minLat = Math.min(minLat, p.lat);
            minLon = Math.min(minLon, p.lon);
            maxLat = Math.max(maxLat, p.lat);
            maxLon = Math.max(maxLon, p.lon);
        }
        return new Rect(minLat, minLon, maxLat, maxLon);
    }

    // ✅ Kiểm tra 2 vùng chữ nhật có giao nhau không
    static boolean intersects(Rect a, Rect b) {
        return !(a.maxLat < b.minLat || a.minLat > b.maxLat
                || a.maxLon < b.minLon || a.minLon > b.maxLon);
    }

    // ✅ Hàm đo thời gian chung
    static <T> T measureTime(String label, Supplier<T> task) {
        long start = System.nanoTime();
        T result = task.get();
        long end = System.nanoTime();
        double elapsedMs = (end - start) / 1_000_000.0;
        System.out.printf("⏰ %s: %.3f ms%n", label, elapsedMs);
        return result;
    }

    public static void main(String[] args) {

        // ====== Khởi tạo dữ liệu lộ trình 10 tài xế ======
         Map<String, List<Point>> drivers = new LinkedHashMap<>();

        // === TX01: chạy hơi chếch tây nam → đông bắc, đi qua khách ===
        drivers.put("TX01", Arrays.asList(
            new Point(20.985, 105.790),
            new Point(20.992, 105.815), // qua khách
            new Point(21.010, 105.860)
        ));

        // === TX02: xa khách ===
        drivers.put("TX02", Arrays.asList(
            new Point(21.015, 105.770),
            new Point(21.040, 105.810)
        ));

        // === TX03: gần như song song TX01, cũng đi qua khách ===
        drivers.put("TX03", Arrays.asList(
            new Point(20.980, 105.805),
            new Point(20.993, 105.816), // qua khách
            new Point(21.020, 105.845)
        ));

        // === TX04: Hồ Tây xa khách ===
        drivers.put("TX04", Arrays.asList(
            new Point(21.030, 105.840),
            new Point(21.060, 105.870)
        ));

        // === TX05: gần Times City, hơi xa khách ===
        drivers.put("TX05", Arrays.asList(
            new Point(21.005, 105.835),
            new Point(21.030, 105.850)
        ));

        // === TX06: Cầu Giấy xa khách ===
        drivers.put("TX06", Arrays.asList(
            new Point(21.045, 105.790),
            new Point(21.070, 105.810)
        ));

        // === TX07: lệch về phía đông, vẫn cắt vùng khách ===
        drivers.put("TX07", Arrays.asList(
            new Point(20.987, 105.812),
            new Point(20.992, 105.832), // gần khách
            new Point(21.015, 105.860)
        ));

        // === TX08: Minh Khai xa khách ===
        drivers.put("TX08", Arrays.asList(
            new Point(21.010, 105.860),
            new Point(21.040, 105.880)
        ));

        // === TX09: Nhật Tân xa khách ===
        drivers.put("TX09", Arrays.asList(
            new Point(21.040, 105.780),
            new Point(21.065, 105.820)
        ));

        // === TX10: Linh Đàm xa khách ===
        drivers.put("TX10", Arrays.asList(
            new Point(20.970, 105.780),
            new Point(21.000, 105.820)
        ));


        // ====== Vị trí khách ======
        Point customer = new Point(20.992, 105.816);

        // ====== Tạo vùng tìm kiếm 1 km quanh khách ======
        Rect searchRect = makeSearchRect(customer, 1.0);

        // ====== Tìm tài xế có lộ trình giao vùng khách ======
        List<String> nearbyDrivers = measureTime("Tìm lộ trình giao vùng 1 km quanh khách", () -> {
            List<String> found = new ArrayList<>();
            for (var entry : drivers.entrySet()) {
                String name = entry.getKey();
                List<Point> route = entry.getValue();
                Rect routeRect = boundingRect(route);
                if (intersects(routeRect, searchRect)) {
                    found.add(name);
                }
            }
            return found;
        });

        // ====== Kết quả ======
        System.out.println("\n=== KẾT QUẢ ===");
        System.out.printf("📍 Vị trí khách: (%.6f, %.6f)%n", customer.lat, customer.lon);
        System.out.printf("🟧 Vùng tìm kiếm: [%.6f, %.6f] → [%.6f, %.6f]%n",
                searchRect.minLat, searchRect.minLon, searchRect.maxLat, searchRect.maxLon);

        if (nearbyDrivers.isEmpty()) {
            System.out.println("❌ Không có tài xế nào giao vùng tìm kiếm.");
        } else {
            System.out.println("✅ Tài xế có lộ trình giao vùng 1 km quanh khách:");
            for (String name : nearbyDrivers) {
                System.out.println("   - " + name);
            }
        }
    }
}
