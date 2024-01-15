#include "scene.hpp"


// Test scene to check correct compilation and run
//  The code should display a rotating cube and a textured quad on top of a circular ground.

using namespace cgp;


void scene_structure::initialize() {
    // Set the behavior of the camera and its initial position
    camera_control.initialize(inputs, window);
    camera_control.set_rotation_axis_z();
    camera_control.look_at({3.0f, 2.0f, 2.0f}, {0, 0, 0}, {0, 0, 1});

    // Create a visual frame representing the coordinate system
    global_frame.initialize_data_on_gpu(mesh_primitive_frame());


    // Create a cube as a mesh
    mesh const cube_mesh = mesh_primitive_cube(/*center*/{0, 0, 0}, /*edge length*/ 1.0f);
    // Initialize a mesh drawable from a mesh structure
    //   - mesh : store buffer of data (vertices, indices, etc) on the CPU. The mesh structure is convenient to manipulate in the C++ code but cannot be displayed (data is not on GPU).
    //   - mesh_drawable : store VBO associated to elements on the GPU + associated uniform parameters. A mesh_drawable can be displayed using the function draw(mesh_drawable, environment). It only stores the indices of the buffers on the GPU - the buffer of data cannot be directly accessed in the C++ code through a mesh_drawable.
    //   Note: a mesh_drawable can be created from a mesh structure in calling [mesh_drawable_name].initialize_data_on_gpu([mesh_name])
    cube.initialize_data_on_gpu(
            cube_mesh);  // note: cube is a mesh_drawable declared as a variable of the class scene_structure
    cube.material.color = {1, 1,
                           0};         // set the color of the cube (R,G,B) - sent as uniform parameter to the shader when display is called

    // Create the ground plane
    ground.initialize_data_on_gpu(mesh_primitive_quadrangle({-2, -2, -1}, {2, -2, -1}, {2, 2, -1}, {-2, 2, -1}));

    // Create the cylinder
    cylinder.initialize_data_on_gpu(
            mesh_primitive_cylinder(/*radius*/ 0.2f, /*first extremity*/{0, -1, 0}, /*second extremity*/{0, 1, 0}));
    cylinder.material.color = {0.8f, 0.8f, 1};

    // Create the sphere
    sphere.initialize_data_on_gpu(mesh_primitive_sphere(0.3f, {0, 0, 0}, 40, 20));
    sphere.material.color = {1.0f, 0.0f, 1.0f};

    // Create a parametric curve
    // **************************************** //
    numarray<vec3> curve_positions;   // the basic structure of a curve is a vector of vec3
    int const N_curve = 150;          // Number of samples of the curve
    for (int k = 0; k < N_curve; ++k) {
        const float u = k / (N_curve - 1.0f); // u \in [0,1]

        // curve oscillating as a cosine
        const float x = 0;

        const float y = 4.0f * (u - 0.5f);
        const float z = 0.1f * std::cos(u * 16 * 3.14f);

        curve_positions.push_back({x, y, z});
    }
    // send data to GPU and store it into a curve_drawable structure
    curve.initialize_data_on_gpu(curve_positions);
    curve.color = {0, 1, 0};

    int N = 100;
    shape = mesh_primitive_grid({0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}, N, N);
    initial_position = shape.position;
    shape_visual.initialize_data_on_gpu(shape);
    // Reset the color of the shape to white (only the texture image will be seen)
    shape_visual.material.color = {1, 1, 1};

// Load the image and associate the texture id to the structure
    shape_visual.texture.load_and_initialize_texture_2d_on_gpu(project::path + "assets/squirrel.jpg");


}


void scene_structure::display_frame() {
    // Set the light to the current position of the camera
    environment.light = camera_control.camera_model.position();

    if (gui.display_frame)
        draw(global_frame, environment);


    float const time = timer.t;
    timer.update();


    // the general syntax to display a mesh is:
    //   draw(mesh_drawable_name, environment);
    //     Note: scene is used to set the uniform parameters associated to the camera, light, etc. to the shader
    draw(ground, environment);

    // Display cylinder
    // ********************************************* //

    // Cylinder rotated around the axis (1,0,0), by an angle = time/2
    vec3 const axis_of_rotation = {1, 1, 0};
    float const angle_of_rotation = 10.0f * time;
    rotation_transform const rotation_cylinder = rotation_transform::from_axis_angle(axis_of_rotation,
                                                                                     angle_of_rotation);

    // Set translation and rotation parameters (these parameters are send and used in shaders through uniform variables when draw is called)
    cylinder.model.rotation = rotation_cylinder;
    cylinder.model.translation = {1.5f, 0, 0};
    draw(cylinder, environment); // Display of the cylinder

    // Meshes can also be displayed as wireframe using the specific draw_wireframe call
    draw_wireframe(cylinder, environment, /*color of the wireframe*/{1, 0.3f, 0.3f});


    // Display cube
    // ********************************************* //

    cube.model.rotation = rotation_transform::from_axis_angle({0, 0, 1}, std::sin(3 * time));
    cube.model.translation = {-1, 0, 0};
    draw(cube, environment);

    // draw sphere
    sphere.model.translation = {0, 2.0f * std::cos(4 * time), 2.0f * std::sin(4 * time)};
    sphere.material.color = vec3(1 + std::cos(time), 1 + std::sin(time), 2.0) / 2.0f;
    draw(sphere, environment);
    if (gui.display_wireframe)
        draw_wireframe(sphere, environment, {1, 1, 0});

    curve.model.translation = {1.9f, 0, 0};
    curve.model.rotation = rotation_transform::from_axis_angle({0, 1, 0}, time);
    draw(curve, environment);

    evolve_shape(time);
    shape_visual.vbo_position.update(shape.position);

    // Recompute normals on the CPU (given the position and the connectivity currently in the mesh structure)
    shape.normal_update();
    // Send updated normals on the GPU
    shape_visual.vbo_normal.update(shape.normal);


    draw(shape_visual, environment);
    if (gui.display_wireframe)
        draw_wireframe(shape_visual, environment, {0, 0, 0});
}

void scene_structure::evolve_shape(float time) {
    size_t const N = initial_position.size();
    for (size_t k = 0; k < N; ++k) {
        vec3 const &p0 = initial_position[k];
        vec3 &p = shape.position[k];
        float const dz = 0.3f * noise_perlin({p0.x + 0.2f * time, p0.y, 0}, 2) +
                         0.015f * noise_perlin({4 * p0.x, 4 * p0.y, time}, 2);
        p = p0 + vec3(0, 0, dz);
    }
}

void scene_structure::display_gui() {
    ImGui::Checkbox("Frame", &gui.display_frame);
    ImGui::SliderFloat("Time Scale", &timer.scale, 0.0f, 2.0f, "%.1f");
    ImGui::Checkbox("Wireframe", &gui.display_wireframe);
}

void scene_structure::mouse_move_event() {
    if (!inputs.keyboard.shift)
        camera_control.action_mouse_move(environment.camera_view);
}

void scene_structure::mouse_click_event() {
    camera_control.action_mouse_click(environment.camera_view);
}

void scene_structure::keyboard_event() {
    camera_control.action_keyboard(environment.camera_view);
}

void scene_structure::idle_frame() {
    camera_control.idle_frame(environment.camera_view);
}

