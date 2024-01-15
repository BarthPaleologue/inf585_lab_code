#include "deformers.hpp"



using namespace cgp;

void apply_deformation(mesh& shape, numarray<vec3> const& position_before_deformation, vec2 const& translate_screen, vec3 const& picked_position, vec3 const& picked_normal, rotation_transform const& camera_orientation, deformer_parameters_structure const& deformer_parameters)
{

	/** 
		shape:  The position of shape are the one to be deformed
		position_before_deformation:  Reference position of the mesh before applying this deformation
		translate_screen:   Input translation of the user in the 2D-screen coordinates - tr must be converted into a 3D transformation applied to the positions of shape
		picked_position:    3D Position of the picked vertex
		picked_normal:      Normal of the surface at the picked vertex position
		camera_orientation: Current camera orientation - allows to convert the 2D-screen coordinates into 3D coordinates
	*/


	float const r = deformer_parameters.falloff; // radius of influence of the deformation
	size_t const N = shape.position.size();
	for (size_t k = 0; k < N; ++k)
	{
		vec3& p_shape = shape.position[k];                             // position to deform
		vec3 const& p_shape_original = position_before_deformation[k]; // reference position before deformation
		vec3 const& p_clicked = picked_position;                       // 3D position of picked position
		vec3 const& n_clicked = picked_normal;                         // normal of the surface (before deformation) at the picked position

		float const dist = norm(p_clicked - p_shape_original);         // distance between the picked position and the vertex before deformation

		// TO DO: Implement the deformation models
		// **************************************************************** //
		// ...
		if (deformer_parameters.type == deform_translate) // Case of translation
		{
			// Hint: You can convert the 2D translation in screen space into a 3D translation in the view plane in multiplying 
			//       camera_orientation * vec3(translate_screen, 0)
			vec3 translation = camera_orientation * vec3(translate_screen, 0.0f);
            if(deformer_parameters.direction == direction_surface_normal) {
                translation = dot(translation, n_clicked) * n_clicked;
            }

			// Fake deformation (linear translation in the screen space) 
			//   the following lines should be modified to get the expected smooth deformation
			if (dist < r) {
                float t = 1.0f - dist / r;
                //float smoothT = 0.5f - 0.5f * cosf(t * M_PI);
                float smoothT = t * t * (3.0f - 2.0f * t);
                p_shape = p_shape_original + smoothT * translation;
            }

		}
		if (deformer_parameters.type == deform_twist)
		{
			vec3 twist_axis = camera_orientation * vec3(0, 0, 1);
            if(deformer_parameters.direction == direction_surface_normal) {
                twist_axis = n_clicked;
            }

            if(dist < r) {
                float angle = translate_screen.x * M_PI * 0.5f * (1.0f - dist / r);
                p_shape = rotation_axis_angle(twist_axis, angle) * (p_shape_original - p_clicked) + p_clicked;
            }
		}
		if (deformer_parameters.type == deform_scale)
		{
			if(dist < r) {
                float t = (1.0f - dist / r);
                float smoothT = t * t * (3.0f - 2.0f * t);
                float scale = 1.0f + translate_screen.x * smoothT;
                p_shape = scale * (p_shape_original - p_clicked) + p_clicked;
            }
		}
        if(deformer_parameters.type == deform_noise) {
            vec3 translation = camera_orientation * vec3(translate_screen, 0.0f);
            if(deformer_parameters.direction == direction_surface_normal) {
                translation = dot(translation, n_clicked) * n_clicked;
            }

            if (dist < r) {
                float t = 1.0f - dist / r;
                //float smoothT = 0.5f - 0.5f * cosf(t * M_PI);
                float smoothT = t * t * (3.0f - 2.0f * t);
                float noise = noise_perlin(p_shape_original, 8, 2.0);
                p_shape = p_shape_original + smoothT * translation * noise / 200.0f;
            }
        }

	}


}