#include "effects.hpp"

using namespace cgp;


// Update a transition_structure before starting a new transition
void effect_transition_start(effect_transition_structure &effect_transition, character_structure &character,
                             std::string const &destination_anim) {
    effect_transition.source_anim = character.current_animation_name;
    effect_transition.timer_source_anim = character.timer;

    effect_transition.destination_anim = destination_anim;
    effect_transition.timer_destination_anim.event_period = character.animated_model.animation.at(
            destination_anim).time_max;
    effect_transition.timer_destination_anim.t_periodic = 0;
    effect_transition.timer_destination_anim.start();

    effect_transition.timer_completion.t = 0;
    effect_transition.timer_completion.scale =
            1.0f / effect_transition.transition_time; // force the timer to go from 0 to 1 during the transition_time
    effect_transition.timer_completion.start();
    effect_transition.active = true;

    character.current_animation_name = destination_anim;
}

// Check if a transition need to be de-activated if it is fully completed
void effect_transition_stop_if_completed(effect_transition_structure &transition, character_structure &character) {
    if (transition.timer_completion.t >= 1.0f) {
        character.set_current_animation(transition.destination_anim);
        character.timer = transition.timer_destination_anim;
        transition.active = false;
    }
}


void effect_transition_compute(effect_transition_structure &effect_transition, character_structure &character) {

    effect_transition_structure &transition = effect_transition;
    animated_model_structure &model = character.animated_model;

    // Compute the skeleton from the source animation
    model.set_skeleton_from_animation(transition.source_anim, transition.timer_source_anim.t_periodic);
    numarray<mat4> joint_source_anim_local = model.skeleton.joint_matrix_local;
    numarray<mat4> joint_source_anim_global = model.skeleton.joint_matrix_local;

    // Compute the skeleton from the destination animation
    model.set_skeleton_from_animation(transition.destination_anim, transition.timer_destination_anim.t_periodic);
    numarray<mat4> joint_destination_anim_local = model.skeleton.joint_matrix_local;
    numarray<mat4> joint_destination_anim_global = model.skeleton.joint_matrix_global;

    float alpha_completion = transition.timer_completion.t;

    /** TO DO : Compute a smooth transition between two animation
     *  The objective is to fill "model.skeleton.joint_matrix_local and global"
     *     with the correct interpolated matrices corresponding to the blending between the source and destination animation
     *
     *
     *
     *  Notes:
     *    - model.skeleton.joint_matrix_local is a numarray<mat4>. Its size corresponds to the number of joints.
     *    - joint_source_anim contains the (local) matrices of the source animation
     *    - joint_destination_anim contains the (local) matrices of the destination animation
     * 	  - alpha_completion contains the ratio of completion.
     *       e.g.: alpha_completion=0 => we have the source anim
     *             alpha_completion=1 => we have the destination anim
     *             alpha_completion=0.5 => we have "0.5 source anim + 0.5 destination anim"
     *    - if model.skeleton.joint_matrix_local is filled, you can update the global matrices using "model.skeleton.update_joint_matrix_local_to_global();"
     *    - if model.skeleton.joint_matrix_global is filled, you can update the global matrices using "model.skeleton.update_joint_matrix_global_to_local();"
     *
     *  Help:
     *    - You can convert a matrix M to its translation t and rotation r (or quaternion q) representation using the following syntax:
     *      affine_rt a = affine_rt::from_matrix(M);
     *      vec3 t = a.translation;
     *      rotation_transform r = a.rotation;
     *      quaternion q = r.get_quaternion();
     *    - An affine_rt structure is a container for a rotation_transform and a translation. Its corresponding 4x4 matrix can be obtained
     *        with {affine_rt}.matrix();
     *
     *
     */

    // Blend the two animations
    for (int i = 0; i < model.skeleton.joint_matrix_local.size(); i++) {
        affine_rt a1 = affine_rt::from_matrix(joint_source_anim_local[i]);
        affine_rt a2 = affine_rt::from_matrix(joint_destination_anim_local[i]);

        vec3 t1 = a1.translation;
        vec3 t2 = a2.translation;

        rotation_transform r1 = a1.rotation;
        rotation_transform r2 = a2.rotation;

        quaternion q1 = r1.get_quaternion();
        quaternion q2 = r2.get_quaternion();

        vec3 t = (1.0f - alpha_completion) * t1 + alpha_completion * t2;

        quaternion q = q1 * (1.0f - alpha_completion) + q2 * alpha_completion;
        q = normalize(q);

        affine_rt a = affine_rt(q, t);

        model.skeleton.joint_matrix_local[i] = a.matrix();
    }

    // Update the global matrices
    model.skeleton.update_joint_matrix_local_to_global();

    // Update the timers
    transition.timer_source_anim.update();
    transition.timer_destination_anim.update();
    transition.timer_completion.update();
    character.timer = transition.timer_destination_anim;
}


// This function implements the change of animation when we press or release the UP key during a walk effect
void effect_walking_keyboard_event(effect_transition_structure &effect_transition, character_structure &character,
                                   cgp::input_devices const &inputs, effect_walking_structure const &effect_walking) {
    // If we just press the key up (or W), start the Walk cycle animation
    if (inputs.keyboard.last_action.is_pressed(GLFW_KEY_UP) || inputs.keyboard.last_action.is_pressed(GLFW_KEY_W)) {
        effect_transition_start(effect_transition, character, effect_walking.walk_anim_name);
    }
    // If we just release the key up (or W), start the Idle cycle animation
    if (inputs.keyboard.last_action.is_released(GLFW_KEY_UP) || inputs.keyboard.last_action.is_released(GLFW_KEY_W)) {
        effect_transition_start(effect_transition, character, effect_walking.idle_anim_name);
    }

}


// This function implements the change of position of the character when a directional key is pressed
// This function is called at every frame when the walk effect is active
void effect_walking(effect_walking_structure &effect_walking, character_structure &character,
                    cgp::input_devices const &inputs, effect_transition_structure const &effect_transition) {

    /** TO DO : An displacement of the character along the direction of the walk.
     *   The character should be able to move straight, and turn, based on the user command.
     *
     *   Inputs:
     *      - The effect_walking structure storing the current position and angle of the character (root joint position and angle)
     *      - The skeleton to be modified by the walk effect
     *      - The current state of the keyboard input (e.g. inputs.keyboard.is_pressed(...))
     *   Output:
     * 		- A modified effect_walking structure and skeleton corresponding to the expected displacement
     *
     *
     *   Help:
     *      - The following syntax allows to check if a specific key is currently pressed:
     *        if( inputs.keyboard.is_pressed({KeyName}) ) {...} , with
     * 		   {KeyName} being: GLFW_KEY_UP, GLFW_KEY_RIGHT, GLFW_KEY_LEFT, GLFW_KEY_DOWN for the arrows
     *                      or: GLFW_KEY_W, GLFW_KEY_A, GLFW_KEY_D, GLFW_KEY_S using for instance the letters
     *      - Given a mat4 structure representing an affine transformation, you can apply the following block operations
     *        - {mat4}.get/set_block_translation({vec3}) : get/set the vec3 block corresponding to the translation part
     *        - {mat4}.get/set_block_linear({mat3}) : get/set the mat3 block corresponding to the linear part
     *      - A rotation with a given axis and angle can be created with the syntax
     *         rotation_transform r = rotation_axis_angle({vec3 axis}, {float angle});
     *         The associated 3x3 matrix can be obtained using r.matrix();
     *      - An internal timer is available in effect_walking.timer
     *        The timer can be updated via effect_walking.timer.update() and the elapsed time since the last update is returned.
     *        This can be used to enforce a constant speed independently of the frame rate
     *
     *
     * 	 Algorithm:
     * 		If(press UP)
     *        Advance straight and update effect_walking.root_position
     *      If(press Left/Right)
     *        Rotate and update effect_walking.root_angle
     *      Update the root joint of the skeleton
     *
     *
     */

    // Update the timer
    effect_walking.timer.update();

    // Get the elapsed time since the last update
    float dt = effect_walking.timer.t;

    // Get the current position and angle of the character
    vec3 &root_position = effect_walking.root_position;
    float &root_angle = effect_walking.root_angle;

    vec3 forward = {0, 0, 1};

    // rotate the forward vector according to the current angle
    rotation_transform r = rotation_axis_angle({0, 1, 0}, root_angle);
    forward = r.matrix() * forward;

    // Compute the displacement of the character
    vec3 displacement = {0, 0, 0};
    if (inputs.keyboard.is_pressed(GLFW_KEY_UP) || inputs.keyboard.is_pressed(GLFW_KEY_W)) {
        displacement += forward * dt * 0.0055;
    }

    float dtheta = 0;

    if(inputs.keyboard.is_pressed(GLFW_KEY_LEFT) || inputs.keyboard.is_pressed(GLFW_KEY_A)) {
        dtheta = dt * 0.01;
    }

    if(inputs.keyboard.is_pressed(GLFW_KEY_RIGHT) || inputs.keyboard.is_pressed(GLFW_KEY_D)) {
        dtheta = -dt * 0.01;
    }

    // update skeleton (for each joint, apply the rotation)
    mat4 &joint_matrix_local = character.animated_model.skeleton.joint_matrix_local[0];
    joint_matrix_local.apply_transform_to_block_linear(rotation_axis_angle({0, 1, 0}, root_angle).matrix());

    vec3 root_initial_position = vec3(0.0f, joint_matrix_local.get_block_translation().y, 0.0f);
    joint_matrix_local.set_block_translation(root_position + root_initial_position);

    // update matrices
    character.animated_model.skeleton.update_joint_matrix_local_to_global();

    effect_walking.root_position += displacement;
    effect_walking.root_angle += dtheta;
}


// Rotate the joint of the head to look toward the given objective position (if possible)
void effect_rotate_head_toward_objective_position(skeleton_structure &skeleton, int joint_head_index,
                                                  vec3 const &objective_position) {
    /** TO DO : Change the direction of the joint of the head such that it points toward the objective position given in parameters
     *
     *  Notes:
     *    - the index of the joint corresponding to the head is supposed to be given in parameter in the variable joint_head_index
     *      The associated matrix is accessible via "skeleton.joint_matrix_local/global[joint_head_index]"
     *    - Given a mat4 structure representing an affine transformation, you can apply the following block operations
     *      - {mat4}.get_block_linear() : returns the mat3 corresponding to the linear part (typically the rotation)
     *      - {mat4}.get_block_translation() : returns the vec3 corresponding to the translation part
     *      - {mat4}.apply_transform_to_block_linear({mat3}) : Applies the mat3 matrix to the linear part of the mat4 matrix (doesn't modify the translation part)
     *    - The transpose of a matrix M (mat3 or mat4) can be obtained with transpose(M), or transpose(M)
     *    - A rotation with a given axis and angle can be created with the syntax
     *      rotation_transform r = rotation_axis_angle({vec3 axis}, {float angle});
     *      The associated 3x3 matrix can be obtained using r.matrix();
     *
     *  skeleton.joint_matrix_local/global[joint_head_index] = ... // something that depends on the objective_position
     *
     */

    mat4 &joint_head_matrix_global = skeleton.joint_matrix_global[15]; // The given index was wrong. So I hardcoded it here.

    vec3 current_forward = {0, 0, 1};
    current_forward = joint_head_matrix_global.get_block_linear() * current_forward;

    vec3 target_forward = normalize(objective_position - joint_head_matrix_global.get_block_translation());

    vec3 rotation_axis = cross(current_forward, target_forward);

    float rotation_angle = acos(dot(current_forward, target_forward));
    rotation_angle = std::min(rotation_angle, 0.9f);
    rotation_angle = std::max(rotation_angle, -0.9f);

    rotation_transform r = rotation_axis_angle(rotation_axis, rotation_angle);

    joint_head_matrix_global.apply_transform_to_block_linear(r.matrix());
}


// Update the current position of the IK constraint to match at its start the current joint position
void effect_ik_start(effect_ik_structure &effect_ik, skeleton_structure &skeleton, int joint_end) {
    // get current position of the joint at index joint_end
    vec3 current_position = skeleton.joint_matrix_global[joint_end].get_block_translation();
    // set this position as the default objective for the IK
    effect_ik.target_position = current_position;
    effect_ik.target_offset = {0, 0, 0};

    for(int i = 0; i < skeleton.joint_matrix_global.size(); i++) {
        int parent = skeleton.parent_index[i];
        effect_ik.parent_index[i] = parent;

        if(parent != -1) {
            effect_ik.child_index[parent] = i;
        }
    }

    int current_joint = joint_end;
    while(current_joint != effect_ik.joint_root_ik) {
        vec3 p1 = skeleton.joint_matrix_global[current_joint].get_block_translation();
        vec3 p2 = skeleton.joint_matrix_global[effect_ik.parent_index[current_joint]].get_block_translation();

        effect_ik.joint_distance_parent[current_joint] = norm(p1 - p2);

        current_joint = effect_ik.parent_index[current_joint];
    }

    current_joint = effect_ik.joint_root_ik;
    while(current_joint != joint_end) {
        vec3 p1 = skeleton.joint_matrix_global[current_joint].get_block_translation();
        vec3 p2 = skeleton.joint_matrix_global[effect_ik.child_index[current_joint]].get_block_translation();

        effect_ik.joint_distance_child[current_joint] = norm(p1 - p2);

        current_joint = effect_ik.child_index[current_joint];
    }
}


void effect_ik_compute(effect_ik_structure const &effect_ik, skeleton_structure &skeleton) {
    /** TO DO : Implement an Inverse Kinematic on the skeleton
     *
     *   Inputs:
     *      - The skeleton to be modified by the IK
     *      - The index of the joint constrained by the IK : effect_ik.joint_end
     *      - The index of the joint at the beginning of the skeleton chain : effect_ik.joint_start
     *        (effect_ik.joint_start should be a parent of effect_ik.joint_end, otherwise the start of the chain can be considered to be the skeleton root)
     *      - The position associated to constrained joint:
     *          p_constrained = effect_ik.objective_position + objective_offset.effect_ik
     *        Rem. this position is split between an initial position and an offset that can be manipulated via the GUI.
     *   Output:
     * 		- A modified joint position of the skeleton satisfying at best the constraint.
     *
     *   Rem.
     *      - effect_ik.joint_start should be a parent of effect_ik.joint_end such that the set of joints between [joint_start - joint_end] is a kinematic chain, subset of the skeleton hierarchy. If it is not the case, we may assume that the start of the chain is at the skeleton root.
     *
     *   Help:
     *      - A rotation R transforming a vector v1 to a vector v2 (with ||v1||=||v2||=1) can be computed using
     *        rotation_transform R = rotation_transform::from_vector_transform(v1, v2);
     *      - Considering a mat4 representing an affine transform:
     *         - {mat4}.get_block_translation() gives the vec3 corresponding to the translation part
     *         - {mat4}.set_block_translation({vec3}) set the block to the mat4 corresponding to the translation with the vec3 passed as argument
     *         - {mat4}.apply_transform_to_block_linear({mat3}) : Applies the mat3 matrix to the linear part of the mat4 matrix (doesn't modify the translation part)
     *         - mat4 M_inv = {mat4}.inverse_assuming_rigid_transform() : gives the 4x4 matrix corresponding to the inverse of the {mat4}
     *      - mat4().set_identity() generate a 4x4 identity matrix
     *
     */

    int joint_end = effect_ik.joint_target;
    int joint_start = effect_ik.joint_root_ik;

    vec3 p_constrained = effect_ik.target_position + effect_ik.target_offset;

    // FabriK algorithm
    // 1. Compute the distance between the constrained joint and the objective position
    vec3 p = skeleton.joint_matrix_global[joint_end].get_block_translation();
    vec3 delta = p_constrained - p;

    // 2. Compute the length of the kinematic chain
    float length = 0;
    int current_joint = joint_end;

    while(current_joint != joint_start) {
        vec3 p1 = skeleton.joint_matrix_global[current_joint].get_block_translation();
        vec3 p2 = skeleton.joint_matrix_global[skeleton.parent_index[current_joint]].get_block_translation();
        length += norm(p1 - p2);
        current_joint = skeleton.parent_index[current_joint];
    }

    // 3. If the distance is greater than the length of the chain, we cannot reach the objective position
    if(norm(delta) > length) {
        return;
    }

    // for a joint, find its parent
    std::map<int, int> parent_index = effect_ik.parent_index;
    // for a joint, find its child in the kinematic chain
    std::map<int, int> child_index = effect_ik.child_index;

    std::map<int, float> joint_distance_parent = effect_ik.joint_distance_parent;
    std::map<int, float> joint_distance_child = effect_ik.joint_distance_child;

    // set position of the constrained joint
    skeleton.joint_matrix_global[joint_end].set_block_translation(p_constrained);

    for(int i = 0; i < 10; i++) {
        // back propagation
        current_joint = joint_end;

        while (current_joint != joint_start) {
            // skip the root joint
            if(current_joint == joint_end) {
                current_joint = parent_index[current_joint];
                continue;
            }

            // get the position of the current joint
            vec3 p1 = skeleton.joint_matrix_global[current_joint].get_block_translation();

            // get the position of the parent joint
            vec3 p2 = skeleton.joint_matrix_global[parent_index[current_joint]].get_block_translation();

            // compute the direction of the joint
            vec3 direction = normalize(p2 - p1);

            // the new position of the joint is the position of the parent joint + the length of the joint * the direction
            vec3 new_position = p2 - joint_distance_parent[current_joint] * direction;

            skeleton.joint_matrix_global[current_joint].set_block_translation(new_position);

            // update the current joint
            current_joint = parent_index[current_joint];

            // update matrices
            skeleton.update_joint_matrix_global_to_local();
        }

        // forward propagation
        current_joint = joint_start;
        while(current_joint != joint_end) {
            // skip the root joint
            if(current_joint == joint_start) {
                current_joint = child_index[current_joint];
                continue;
            }

            // get the position of the current joint
            vec3 p1 = skeleton.joint_matrix_global[current_joint].get_block_translation();
            // get the position of the child joint
            vec3 p2 = skeleton.joint_matrix_global[child_index[current_joint]].get_block_translation();

            // compute the direction of the joint
            vec3 direction = normalize(p2 - p1);

            // the new position of the joint is the position of the child joint - the length of the joint * the direction
            vec3 new_position = p2 - joint_distance_child[current_joint] * direction;

            skeleton.joint_matrix_global[current_joint].set_block_translation(new_position);

            // update the current joint
            current_joint = child_index[current_joint];

            // update matrices
            skeleton.update_joint_matrix_global_to_local();
        }
    }

    // move children of end joint by delta
    std::function<void(int)> moveJoint = [&](int joint) {
        vec3 p = skeleton.joint_matrix_global[joint].get_block_translation();
        skeleton.joint_matrix_global[joint].set_block_translation(p + delta);

        for(int child : skeleton.child(joint)) {
            moveJoint(child);
        }
    };

    for(int child : skeleton.child(joint_end)) {
        moveJoint(child);
    }
}


