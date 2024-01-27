#include "animated_model.hpp"
#include "DualQuaternion.h"


using namespace cgp;


void animated_model_structure::skinning_lbs() {
    // ************************************************************** //
    // TO DO: Compute the Linear Blend Skinning (LBS) deformation
    // ...
    // ************************************************************** //
    //
    // Help:
    //     - The function should update the values of rigged_mesh.mesh_deformed.position based on the skeleton and bind pose (rigged_mesh.mesh_bind_pose.position)
    //     - Once the computation is working on the position, you may also update the normals
    //     - The skinning weights are available via: rigged_mesh.skinning_weight
    //       They are stored per vertex and per joint: float weight_ij = rigged_mesh.skinning_weight[vertex_i][joint_j];
    //       
    //     - Given a mat4 M representing a rigid transformation (with rotation and translation only), you can compute efficiently its inverse using the syntax
    //        mat4 M_inversed = M.inverse_assuming_rigid_transform();
    //     - Consider a mat4 M representing a projective (or affine, or rigid) transformation, and a vec3 p. We call q the transformation of p by M.
    //         - If p represents a 3D point, then q can be expressed as
    //             vec3 q = vec3( M * vec4(p,1.0f) ); or similarily vec3 q = M.transform_position(p);
    //         - If p represents a 3D vector, then q can be expressed as
    //             vec3 q = vec3( M * vec4(p,0.0f) ); or similarily vec3 q = M.transform_vector(p);
    //  

    // Example of looping over the positions of the mesh
    int N_vertex = rigged_mesh.mesh_bind_pose.position.size();
    for (int k_vertex = 0; k_vertex < N_vertex; ++k_vertex) {
        vec3 const &position_in_bind_pose = rigged_mesh.mesh_bind_pose.position[k_vertex]; // The "initial/bind pose" position p0
        vec3 const &normal_in_bind_pose = rigged_mesh.mesh_bind_pose.normal[k_vertex];     // The "initial/bind pose" normal n0
        vec3 &position_to_be_deformed = rigged_mesh.mesh_deformed.position[k_vertex];      // The position to be deformed by LBS
        vec3 &normal_to_be_deformed = rigged_mesh.mesh_deformed.normal[k_vertex];         // The normal to be deformed by LBS

        position_to_be_deformed = vec3();
        normal_to_be_deformed = vec3();

        int nb_joints = rigged_mesh.skinning_weight[k_vertex].size();
        for(int k_joint = 0; k_joint < nb_joints; k_joint++) {
            float weigh_ij = rigged_mesh.skinning_weight[k_vertex][k_joint];

            auto M = skeleton.joint_matrix_global[k_joint];
            auto M0 = skeleton.joint_matrix_global_bind_pose[k_joint];
            auto invM0 = M0.inverse_assuming_rigid_transform();
            auto T = M * invM0;

            position_to_be_deformed += weigh_ij * T.transform_position(position_in_bind_pose);
            normal_to_be_deformed += weigh_ij * T.transform_vector(normal_in_bind_pose);
        }

        // Do some computation ...
        //position_to_be_deformed = position_in_bind_pose;   // to be changed
        //normal_to_be_deformed = normal_in_bind_pose; // to be changed
    }

}

void animated_model_structure::skinning_dqs() {
    // ************************************************************** //
    // TO DO: Compute Dual Quaternion Skinning (DQS) deformation
    // ...
    // ************************************************************** //
    //
    // Help:
    //     - Given a mat4 representing a rigid transformation, the following syntax allows to access the rotation and translation part:
    //         affine_rt a = affine_rt::from_matrix({mat4});
    //         rotation_transform rot = a.rotation
    //         vec3 translation = a.translation
    //     - The quaternion of a rotation_transform can be accessed via {rotation_transform}.get_quaternion();
    //     - The structure quaternion is a specialized type derived from a vec4. You can access to its .x .y .z .w component similarily to a vec4.
    //

    int N_vertex = rigged_mesh.mesh_bind_pose.position.size();
    for(int k_vertex = 0; k_vertex < N_vertex; k_vertex++) {
        vec3 const &position_in_bind_pose = rigged_mesh.mesh_bind_pose.position[k_vertex]; // The "initial/bind pose" position p0
        vec3 const &normal_in_bind_pose = rigged_mesh.mesh_bind_pose.normal[k_vertex];     // The "initial/bind pose" normal n0
        vec3 &position_to_be_deformed = rigged_mesh.mesh_deformed.position[k_vertex];      // The position to be deformed by LBS
        vec3 &normal_to_be_deformed = rigged_mesh.mesh_deformed.normal[k_vertex];         // The normal to be deformed by LBS

        auto dualQuaternionSum = DualQuaternion::FromRotationTranslation(quaternion(0,0,0,0), vec3(0,0,0));

        int nb_joints = rigged_mesh.skinning_weight[k_vertex].size();
        for(int k_joint = 0; k_joint < nb_joints; k_joint++) {
            auto weigh_ij = rigged_mesh.skinning_weight[k_vertex][k_joint];

            auto M = skeleton.joint_matrix_global[k_joint];
            auto rotationQuaternion = affine_rt::from_matrix(M).rotation.get_quaternion();
            auto translation = affine_rt::from_matrix(M).translation;

            auto dualQuaternion = DualQuaternion::FromRotationTranslation(rotationQuaternion, translation);

            dualQuaternionSum = dualQuaternionSum + dualQuaternion * weigh_ij;
        }

        auto transformation = dualQuaternionSum.getMatrix();

        position_to_be_deformed = transformation.transform_position(position_in_bind_pose);
        normal_to_be_deformed = transformation.transform_vector(normal_in_bind_pose);
    }
}
