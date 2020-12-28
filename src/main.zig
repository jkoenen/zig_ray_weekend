// Ray Tracing in One Weekend:
// https://raytracing.github.io/books/RayTracingInOneWeekend.html#overview
// lets test zig..

const std = @import("std");
const math = std.math;

const vec3 = @import("vec3.zig");
const raylib = @import("ray.zig");
const Color = vec3.Color;
const Vec3 = vec3.Vec3;
const Point = vec3.Point;
const Ray = raylib.Ray;

fn createNormalizedPosition(x: i32, y: i32, width: i32, height: i32, z: f32) Vec3 {
    // is there a nicer / simpler way for this?
    const xf = @intToFloat(f32, x) / @intToFloat(f32, width - 1);
    const yf = @intToFloat(f32, y) / @intToFloat(f32, height - 1);
    const zf = z;

    return .{ .x = xf, .y = yf, .z = zf };
}

fn skyColor(r: Ray) Color {
    const unit_direction = vec3.normalize(r.direction);
    const t = 0.5 * (unit_direction.y + 1.0);
    const topColor = Color.new(1, 1, 1);
    return vec3.lerp(Color.new(0.5, 0.7, 1), Color.new(1, 1, 1), t);
}

fn hitSphere(center: Point, radius: f32, ray: Ray) f32 {
    const oc = vec3.sub(ray.origin, center);
    const a = vec3.dot(ray.direction, ray.direction);
    const b = 2 * vec3.dot(oc, ray.direction);
    const c = vec3.dot(oc, oc) - radius * radius;
    const discriminant = b * b - 4 * a * c;
    if (discriminant < 0) {
        return -1;
    } else {
        return (-b - math.sqrt(discriminant)) / (2 * a);
    }
}

fn rayColor(r: Ray) Color {
    const sphere_center = Point.new(0, 0, -1);
    const sphere_radius = 0.5;
    const t = hitSphere(sphere_center, sphere_radius, r);
    if (t > 0) {
        const normal = vec3.normalize(vec3.sub(r.at(t), sphere_center));
        return Color.new(0.5 * (normal.x + 1), 0.5 * (normal.y + 1), 0.5 * (normal.z + 1));
    }
    return skyColor(r);
}

fn writePpmColor(out: anytype, color: Color) !void {
    const ir = @floatToInt(i32, 255.0 * color.x + 0.5);
    const ig = @floatToInt(i32, 255.0 * color.y + 0.5);
    const ib = @floatToInt(i32, 255.0 * color.z + 0.5);

    try std.fmt.format(out, "{} {} {}\n", .{ ir, ig, ib });
}

pub fn main() anyerror!void {
    const stdout = std.io.getStdOut().writer();

    // first step: output a ppm:
    const aspect_ratio: f32 = 16.0 / 9.0;
    const image_width: i32 = 256;
    const image_height: i32 = 256 / aspect_ratio;

    const viewport_height: f32 = 2.0;
    const viewport_width: f32 = aspect_ratio * viewport_height;
    const focal_length: f32 = 1;

    const origin = Vec3.new(0, 0, 0);
    const horizontal = Vec3.new(viewport_width, 0, 0);
    const vertical = Vec3.new(0, viewport_height, 0);
    const lower_left_corner = Vec3.new(origin.x - horizontal.x / 2 - vertical.x / 2, origin.y - horizontal.y / 2 - vertical.y / 2, origin.z - focal_length);

    const ppm_file = try std.fs.cwd().createFile("test.ppm", .{});
    defer ppm_file.close();

    const ppm_writer = ppm_file.writer();

    try ppm_writer.print("P3\n{} {}\n255\n", .{ image_width, image_height });

    // hmm.. how do for-loops work in zig? -> using while loops for now ;)
    var y: i32 = image_height - 1;
    while (y >= 0) {
        try stdout.print("\rScanlines remaining: {}", .{y});
        var x: i32 = 0;
        while (x < image_width) {
            const u = @intToFloat(f32, x) / @intToFloat(f32, image_width - 1);
            const v = @intToFloat(f32, y) / @intToFloat(f32, image_height - 1);

            const ray = Ray.new(origin, Vec3.new(lower_left_corner.x + u * horizontal.x + v * vertical.x - origin.x, lower_left_corner.y + u * horizontal.y + v * vertical.y - origin.y, lower_left_corner.z + u * horizontal.z + v * vertical.z - origin.z));

            const pixel_color = rayColor(ray);
            try writePpmColor(ppm_writer, pixel_color);

            x += 1;
        }
        y -= 1;
    }

    try stdout.print("\nDone\n", .{});
}
