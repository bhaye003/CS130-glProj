#include "driver_state.h"
#include <cstring>

driver_state::driver_state()
{
}

driver_state::~driver_state()
{
    delete [] image_color;
    delete [] image_depth;
}

// This function should allocate and initialize the arrays that store color and
// depth.  This is not done during the constructor since the width and height
// are not known when this class is constructed.
void initialize_render(driver_state& state, int width, int height)
{
    	state.image_width=width;
    	state.image_height=height;
    	state.image_color=0;
   	state.image_depth=0;

    	int length = width * height;
	state.image_color = new pixel[length];
	
	for(int i = 0; i < length; i++){
		state.image_color[i] = make_pixel(0,0,0);
	}
	state.image_depth = new float[length];
	for(int i = 0; i < length; i++){
		state.image_depth[i] = 1;
	}
}

// This function will be called to render the data that has been stored in this class.
// Valid values of type are:
//   render_type::triangle - Each group of three vertices corresponds to a triangle.
//   render_type::indexed -  Each group of three indices in index_data corresponds
//                           to a triangle.  These numbers are indices into vertex_data.
//   render_type::fan -      The vertices are to be interpreted as a triangle fan.
//   render_type::strip -    The vertices are to be interpreted as a triangle strip.
void render(driver_state& state, render_type type)
{
switch(type){
	case render_type::triangle: {
		const data_geometry* arr[3];
		data_geometry datageo[3];
		data_vertex dataver[3];
		int tri_ver = state.num_vertices / 3;
		int VDindex = 0;
		
		for(int i = 0; i < tri_ver; i++){
			for(int j = 0; j < 3; j++){
				dataver[j].data = &state.vertex_data[VDindex];
				datageo[j].data = dataver[j].data;
				state.vertex_shader(dataver[j], datageo[j], state.uniform_data);
				arr[j] = &datageo[j];
				VDindex += state.floats_per_vertex;
			}
			clip_triangle(state, arr, 0);
		}
		break;
	}
	case render_type::indexed:{
		const data_geometry* arr[3];
		data_geometry datageo[3];
		data_vertex dataver[3];
		
		for(int i = 0; i < 3 * state.num_triangles; i += 3){
			for(int j = 0; j < 3; j++){
				dataver[j].data = &state.vertex_data[state.index_data[i + j] * state.floats_per_vertex];
				datageo[j].data = dataver[j].data;
				state.vertex_shader(dataver[j], datageo[j], state.uniform_data);
				arr[j] = &datageo[j];
			}
			clip_triangle(state, arr, 0);
		}
		break;
	}
	case render_type::fan: {
		const data_geometry* arr[3];
		data_geometry datageo[3];
		data_vertex dataver[3];
	
		for(int i = 0; i < state.num_vertices; i++){
			for(int j = 0; j < 3; j++){
				int index = i + j;
				if(j == 0){
					index = 0;
				}
				dataver[j].data = &state.vertex_data[index * state.floats_per_vertex];
				datageo[j].data = dataver[j].data;
				state.vertex_shader(dataver[j],datageo[j], state.uniform_data);
				arr[j] = &datageo[j];
			}
			clip_triangle(state, arr, 0);
		}
		break;
	}
	case render_type::strip: {
		const data_geometry* arr[3];
		data_geometry datageo[3];
		data_vertex dataver[3];

		for(int i = 0; i < state.num_vertices - 2; i++){
			for(int j = 0; j<3; j++){
				dataver[j].data = &state.vertex_data[(i+j) * state.floats_per_vertex];
				datageo[j].data = dataver[j].data;
				state.vertex_shader(dataver[j], datageo[j], state.uniform_data);
				arr[j] = &datageo[j];
			}
			clip_triangle(state, arr, 0);
		}
		break;
	}
	default:{
		break;
	}
}
}

// This function clips a triangle (defined by the three vertices in the "in" array).
// It will be called recursively, once for each clipping face (face=0, 1, ..., 5) to
// clip against each of the clipping faces in turn.  When face=6, clip_triangle should
// simply pass the call on to rasterize_triangle.
void clip_triangle(driver_state& state, const data_geometry* in[3],int face)
{

    if(face==1)
    {
        rasterize_triangle(state, in);
        return;
    }


    

    vec4 a = in[0]->gl_Position;
    vec4 b = in[1]->gl_Position;
    vec4 c = in[2]->gl_Position;

    const data_geometry *input[3] = {in[0], in[1], in[2]};
    data_geometry input_dataA[3];
    data_geometry input_dataB[3];
    
    float fa;
    float fb0;
    float fb1;
    
    vec4 p0, p1;
 
    if (a[2] < -a[3] && b[2] < -b[3] && c[2] < -c[3]){

        return;
    }
    
    else if (a[2] < -a[3] && b[2] >= -b[3] && c[2] >= -c[3]){

            fb0 = (-b[3] - b[2]) / (a[2] + a[3] - b[3] - b[2]);
            fb1 = (-a[3] - a[2]) / (c[2] + c[3] - a[3] - a[2]);
            p0 = fb0 * a + (1 - fb0) * b;
            p1 = fb1 * c + (1 - fb1) * a;

            input_dataA[0].data = new float[state.floats_per_vertex];
            input_dataA[1] = *in[1];
            input_dataA[2] = *in[2];

            for (int i = 0; i < state.floats_per_vertex; ++i){
                switch (state.interp_rules[i])
                {
                case interp_type::flat:
                    input_dataA[0].data[i] = in[0]->data[i];
                    break;
                case interp_type::smooth:
                    input_dataA[0].data[i] = fb1 * in[2]->data[i] + (1 - fb1) * in[0]->data[i];
                    break;
                case interp_type::noperspective:
                    fa = fb1 * in[2]->gl_Position[3] / (fb1 * in[2]->gl_Position[3] + (1 - fb1) * in[0]->gl_Position[3]);
                    input_dataA[0].data[i] = fa * in[2]->data[i] + (1 - fa) * in[0]->data[i];
                    break;
                default:
                    break;
                }
	    }

            input_dataA[0].gl_Position = p1;
            input[0] = &input_dataA[0];
            input[1] = &input_dataA[1];
            input[2] = &input_dataA[2];
            
            clip_triangle(state, input, face + 1);

            input_dataB[0].data = new float[state.floats_per_vertex];
            input_dataB[2] = *in[2];

            for (int i = 0; i < state.floats_per_vertex; ++i){

                switch (state.interp_rules[i])
                {
                case interp_type::flat:
                    input_dataB[0].data[i] = in[0]->data[i];
                    break;
                case interp_type::smooth:
                    input_dataB[0].data[i] = fb0 * in[0]->data[i] + (1 - fb0) * in[1]->data[i];
                    break;
                case interp_type::noperspective:
                    fa = fb0 * in[0]->gl_Position[3] / (fb0 * in[0]->gl_Position[3] + (1 - fb0) * in[1]->gl_Position[3]);
                    input_dataB[0].data[i] = fa * in[0]->data[i] + (1 - fa) * in[1]->data[i];
                    break;
                default:
                    break;
                }
	    }

            input_dataB[0].gl_Position = p0;
            input[0] = &input_dataB[0];
            input[1] = &input_dataA[1];
	    input[2] = &input_dataA[0];

	}

clip_triangle(state,input,face+1);
}

// Rasterize the triangle defined by the three vertices in the "in" array.  This
// function is responsible for rasterization, interpolation of data to
// fragments, calling the fragment shader, and z-buffering.
void rasterize_triangle(driver_state& state, const data_geometry* in[3])
{
	int iwt = state.image_width;
	int iht = state.image_height;

	float ax = 0;
	float ay = 0;
	float bx = 0;
	float by = 0;
	float cx = 0;
	float cy = 0;
	
	int px = 0;
	int py = 0;
	
	float abcA = 0;
	float pbcA = 0;
	float apcA = 0;
	float abpA = 0;
	
	float alpha = 0;
	float beta = 0;
	float gamma = 0;

	float aper = 0;
	float bper = 0;
	float gper = 0;
	
	ax = (iwt/2)*(in[0]->gl_Position[0] / in[0]->gl_Position[3]) + (iwt/2) - (.5);
	ay = (iht/2)*(in[0]->gl_Position[1] / in[0]->gl_Position[3]) + (iht/2) - (.5);
	
	bx = (iwt/2)*(in[1]->gl_Position[0] / in[1]->gl_Position[3]) + (iwt/2) - (.5);
	by = (iht/2)*(in[1]->gl_Position[1] / in[1]->gl_Position[3]) + (iht/2) - (.5);
	
	cx = (iwt/2)*(in[2]->gl_Position[0] / in[2]->gl_Position[3]) + (iwt/2) - (.5);
	cy = (iht/2)*(in[2]->gl_Position[1] / in[2]->gl_Position[3]) + (iht/2) - (.5);	
	
	
	float minx = std::min(ax, std::min(bx, cx));
	float miny = std::min(ay, std::min(by, cy));
	float maxx = std::max(ax, std::max(bx, cx));
	float maxy = std::max(ay, std::max(by, cy));

	abcA = 0.5 * (((bx*cy)-(cx*by)) - ((ax*cy)-(cx*ay)) + ((ax*by)-(bx*ay)));

	if(minx < 0){
		minx = 0;
	}
	if(maxx > state.image_width){
		maxx = state.image_width;
	}
	if(miny < 0){
		miny = 0;
	}
	if(maxy > state.image_height){
		maxy = state.image_height;
	}

	for(px = minx; px < maxx; px++){
		for(py = miny; py < maxy; py++){
			int index = px + py * state.image_width;
			pbcA = 0.5 * (((bx*cy) - (cx*by)) + ((by - cy) * px) + ((cx - bx )*py));
			apcA = 0.5 * (((cx*ay) - (ax*cy)) + ((cy - ay) * px) + ((ax  - cx )*py));
			abpA = 0.5 * (((ax*by) - (bx*ay)) + ((ay - by) * px) + ((bx - ax )*py));
			
			alpha = pbcA / abcA;
			beta = apcA / abcA;
			gamma = abpA / abcA;
			

			if(alpha >= 0 && beta >= 0 && gamma >= 0){
				data_fragment frag;
				frag.data = new float[MAX_FLOATS_PER_VERTEX];
				data_output o;
				float z_depth = (alpha * (in[0]->gl_Position[2] / in[0]->gl_Position[3])) + (beta * (in[1]->gl_Position[2] / in[1]->gl_Position[3])) + (gamma * (in[2]->gl_Position[2] / in[2]->gl_Position[3]));
				if(state.image_depth[index] > z_depth){
					for(int float_index = 0; float_index < state.floats_per_vertex; float_index++){
						float fl;
						switch(state.interp_rules[float_index]){
							case interp_type::flat:
								frag.data[float_index] = in[0] -> data[float_index];
								break;
							case interp_type::smooth:
								fl = (alpha/in[0]->gl_Position[3]) + (beta/in[1]->gl_Position[3]) + (gamma/in[2]->gl_Position[3]);
								aper = alpha/ (in[0]->gl_Position[3] * fl);
								bper = beta / (in[1]->gl_Position[3] * fl);
								gper = gamma / (in[2]->gl_Position[3] * fl);
								frag.data[float_index] = (aper * in[0]->data[float_index]) + (bper * in[1]->data[float_index]) + (gper * in[2]->data[float_index]);
								break;
							case interp_type::noperspective:
								frag.data[float_index] = alpha*in[0]->data[float_index] + beta*in[1]->data[float_index] + gamma*in[2]->data[float_index];
								break;
							default:
								break;
							}
						}
					state.fragment_shader(frag, o, state.uniform_data);
					state.image_color[index] = make_pixel(o.output_color[0] * 255, o.output_color[1] * 255, o.output_color[2]*255);
					state.image_depth[index] = z_depth;
				}
			}
		}
	}	

    
}
