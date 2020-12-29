// Ray Tracing in One Weekend:
// https://raytracing.github.io/books/RayTracingInOneWeekend.html#overview
// lets test zig..

// open questions:
// - how to switch between debug/release build configs?
//      -> debug is default, -O ReleaseFast switches to optimized build

const std = @import("std");

usingnamespace @import("math.zig");
usingnamespace @import("ray.zig");

fn writePpmColor(out: anytype, color: Color) !void {
    const g_r = std.math.sqrt(color.x);
    const g_g = std.math.sqrt(color.y);
    const g_b = std.math.sqrt(color.z);
    const ir = clamp(@floatToInt(i32, 255.0 * g_r + 0.5), 0, 255);
    const ig = clamp(@floatToInt(i32, 255.0 * g_g + 0.5), 0, 255);
    const ib = clamp(@floatToInt(i32, 255.0 * g_b + 0.5), 0, 255);

    try std.fmt.format(out, "{} {} {}\n", .{ ir, ig, ib });
}

pub fn main() anyerror!void {
    const stdout = std.io.getStdOut().writer();

    // allocator:
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer {
        const leaked = gpa.deinit();
        if (leaked) {
            const result = std.io.getStdErr().writer().print( "Memory leaks detected!\n", .{} );
        }
    }
    var allocator = &gpa.allocator;

    // init random numbers:
    const time_now = std.time.timestamp();
    const seed_i64 = if (time_now >= 0) time_now else -time_now;
    const seed_u64 = @intCast(u64, seed_i64);
    var rng = std.rand.DefaultPrng.init(seed_u64);
    random_init(&rng.random);

    // create the scene/scene:
    const sphere_capacity = 128;
    const material_capacity = 128;
    var scene = try RayScene.init(allocator, sphere_capacity, material_capacity);
    defer scene.deinit();

    const sceneId = 3;
    if (sceneId == 0) {
        const matGreen = scene.addLambertMaterial(Color.new(0.2, 0.6, 0.1));
        const matRed = scene.addLambertMaterial(Color.new(0.7, 0.2, 0.1));
        const matMetal = scene.addMetalMaterial(Color.new(0.7, 0.7, 0.7),0.0);

        scene.addSphere(Point3.new(0, 0, -1), 0.5, matRed);
        scene.addSphere(Point3.new(0, 0.7, -1), 0.3, matMetal);
        scene.addSphere(Point3.new(0, 1.1, -1), 0.2, matRed);
        scene.addSphere(Point3.new(0, -100.5, -1), 100, matGreen);
    } else if (sceneId == 1) {
        const material_ground = scene.addLambertMaterial(Color.new(0.8, 0.8, 0.0));
        const material_center = scene.addLambertMaterial(Color.new(0.1, 0.2, 0.5));
        const material_left = scene.addDielectricMaterial(1.5);
        const material_right = scene.addMetalMaterial(Color.new(0.8, 0.6, 0.2), 0.0);

        scene.addSphere(Point3.new(0.0, -100.5, -1.0), 100.0, material_ground);
        scene.addSphere(Point3.new(0.0, 0.0, -1.0), 0.5, material_center);
        scene.addSphere(Point3.new(-1.0, 0.0, -1.0), 0.5, material_left);
        scene.addSphere(Point3.new(-1.0, 0.0, -1.0), -0.4, material_left);
        scene.addSphere(Point3.new(1.0, 0.0, -1.0), 0.5, material_right);
    } else if (sceneId == 2) {
        const material_left = scene.addLambertMaterial(Color.new(0, 0, 1));
        const material_right = scene.addLambertMaterial(Color.new(1, 0, 0));

        const radius = std.math.cos(@as(f32, std.math.pi / 4.0));
        scene.addSphere(Point3.new(-radius, 0.0, -1.0), radius, material_left);
        scene.addSphere(Point3.new(radius, 0.0, -1.0), radius, material_right);
    } else if (sceneId == 3) {
        const material_ground = scene.addLambertMaterial(Color.new(0.8, 0.8, 0.0));
        const material_center = scene.addLambertMaterial(Color.new(0.1, 0.2, 0.5));
        const material_left = scene.addDielectricMaterial(1.5);
        const material_right = scene.addMetalMaterial(Color.new(0.8, 0.6, 0.2), 0.0);

        scene.addSphere(Point3.new(0.0, -100.5, -1.0), 100.0, material_ground);
        scene.addSphere(Point3.new(0.0, 0.0, -1.0), 0.5, material_center);
        scene.addSphere(Point3.new(-1.0, 0.0, -1.0), 0.5, material_left);
        scene.addSphere(Point3.new(-1.0, 0.0, -1.0), -0.45, material_left);
        scene.addSphere(Point3.new(1.0, 0.0, -1.0), 0.5, material_right);
    }

    const resolution_scale: i32 = 4;
    const samples_per_pixel: i32 = 10;
    const max_depth: i32 = 50;

    const image_width = 1280 / resolution_scale;
    const image_height = 720 / resolution_scale;
    const aspect_ratio = @intToFloat(f32, image_width) / @intToFloat(f32, image_height);

    // camera settings:
    const look_from = Point3.new(3, 3, 2);
    const look_at = Point3.new(0, 0, -1);
    const up = Vector3.new(0, 1, 0);
    const distance_to_focus = sub(look_from, look_at).length();
    const aperture = 2.0;
    const camera = Camera.new(look_from, look_at, up, 20.0, aspect_ratio, aperture, distance_to_focus);

    // open the output file:
    const ppm_file = try std.fs.cwd().createFile("test.ppm", .{});
    defer ppm_file.close();

    const ppm_writer = ppm_file.writer();

    try ppm_writer.print("P3\n{} {}\n255\n", .{ image_width, image_height });

    // hmm.. how do for-loops work in zig? -> using while loops for now ;)
    var y: i32 = image_height - 1;
    while (y >= 0) {
        try stdout.print("\rScanlines remaining: {} ", .{y});
        var x: i32 = 0;
        while (x < image_width) {
            var pixel_color = Color{ .x = 0, .y = 0, .z = 0 };

            var sample_index: i32 = 0;
            while (sample_index < samples_per_pixel) {
                const u = (@intToFloat(f32, x) + rng.random.float(f32)) / @intToFloat(f32, image_width - 1);
                const v = (@intToFloat(f32, y) + rng.random.float(f32)) / @intToFloat(f32, image_height - 1);

                const ray = camera.getRay(u, v);

                const sample_color = scene.rayColor( ray, max_depth);

                pixel_color.add(sample_color);
                sample_index += 1;
            }
            pixel_color.scale(1.0 / @intToFloat(f32, samples_per_pixel));
            try writePpmColor(ppm_writer, pixel_color);

            x += 1;
        }
        y -= 1;
    }

    try stdout.print("\nDone\n", .{});
}
