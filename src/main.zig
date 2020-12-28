// Ray Tracing in One Weekend:
// https://raytracing.github.io/books/RayTracingInOneWeekend.html#overview
// lets test zig..

const std = @import("std");

const vec3 = @import("vec3.zig");
const Color = vec3.Color;
const Vec3 = vec3.Vec3;

fn createNormalizedPosition(x: i32, y: i32, width: i32, height: i32, z: f32) Vec3 {
    // is there a nicer / simpler way for this?
    const xf = @intToFloat(f32, x) / @intToFloat(f32, width - 1);
    const yf = @intToFloat(f32, y) / @intToFloat(f32, height - 1);
    const zf = z;

    return .{ .x = xf, .y = yf, .z = zf };
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
    const image_width: i32 = 256;
    const image_height: i32 = 256;

    var p0: Vec3 = .{ .x = 1, .y = 2, .z = 3 };

    const l: f32 = p0.length();

    try stdout.print("l={}\n", .{l});

    //unreachable;

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
            const pixel_color = createNormalizedPosition(x, y, image_width, image_height, 0.25);
            try writePpmColor(ppm_writer, pixel_color);

            x += 1;
        }
        y -= 1;
    }

    try stdout.print("\nDone\n", .{});
}
