// Ray Tracing in One Weekend:
// https://raytracing.github.io/books/RayTracingInOneWeekend.html#overview
// lets test zig..

const std = @import("std");

const vec3 = @import("vec3.zig");

pub fn main() anyerror!void {
    const stdout = std.io.getStdOut().writer();

    // first step: output a ppm:
    const image_width: i32 = 256;
    const image_height: i32 = 256;

    var p0: vec3.Vec3 = .{ .x = 1, .y = 2, .z = 3 };

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
            // is there a nicer / simpler way for this?
            const r = @intToFloat(f32, x) / @intToFloat(f32, image_width - 1);
            const g = @intToFloat(f32, y) / @intToFloat(f32, image_height - 1);
            const b: f32 = 0.25;

            const ir = @floatToInt(i32, 255.0 * r + 0.5);
            const ig = @floatToInt(i32, 255.0 * g + 0.5);
            const ib = @floatToInt(i32, 255.0 * b + 0.5);

            try ppm_writer.print("{} {} {}\n", .{ ir, ig, ib });

            x += 1;
        }
        y -= 1;
    }

    try stdout.print("\nDone\n", .{});
}
