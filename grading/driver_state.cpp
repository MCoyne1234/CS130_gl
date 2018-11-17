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
    unsigned int area = width*height; //MC

    state.image_width=width;
    state.image_height=height;
    
    state.image_color = new pixel[area];// 0;
    state.image_depth=0;
    //std::cout<<"TODO: allocate and initialize state.image_color and state.image_depth."<<std::endl;
    pixel black = make_pixel(0,0,0);
    for(unsigned int i = 0; i < area ; ++i){
        state.image_color[i] = black;
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
    //std::cout<<"TODO: implement rendering."<<std::endl;
    switch (type){
        case render_type::triangle :
            {
            /*
            for(int verts = 0 ;verts < state.num_vertices/3; ++verts){
                int iter = verts * 9;
                data_geometry dg_arr[3];
                const data_geometry * p_dg_arr[3] = {&dg_arr[0], &dg_arr[1], &dg_arr[2]};
                for(int i = 0; i < 3 ; ++i ){
                    int incr = iter + i * state.floats_per_vertex ;
                    dg_arr[i].data = (state.vertex_data+incr);
                    }
                    rasterize_triangle(state, p_dg_arr);
            }*/
            data_geometry dg_arr[3];
            const data_geometry * p_dg_arr[3] = {&dg_arr[0], &dg_arr[1], &dg_arr[2]};
            int incr = 3 * state.floats_per_vertex;
            for(int i = 0; i < state.num_vertices*state.floats_per_vertex; i += incr){
               // std::cout << "YO " << state.num_vertices << std::endl;
                    dg_arr[0].data = (state.vertex_data+i);
                    dg_arr[1].data = (state.vertex_data+i+state.floats_per_vertex);
                    dg_arr[2].data = (state.vertex_data+i+2*state.floats_per_vertex);
                    rasterize_triangle(state, p_dg_arr);
            }
            }
            break;
        case render_type::indexed :
            {
            std::cout << "indexed" << std::endl;
            break;
            }
        case render_type::fan :
            {
            std::cout << "fan" << std::endl;
            break;
            }
        case render_type::strip :
            {
            std::cout << "strip" << std::endl;
            break;
            }
        default:
            {
            std::cout << "BAD!" << std::endl;
            }
    }
}


// This function clips a triangle (defined by the three vertices in the "in" array).
// It will be called recursively, once for each clipping face (face=0, 1, ..., 5) to
// clip against each of the clipping faces in turn.  When face=6, clip_triangle should
// simply pass the call on to rasterize_triangle.
void clip_triangle(driver_state& state, const data_geometry* in[3],int face)
{
    if(face==6)
    {
        rasterize_triangle(state, in);
        return;
    }
    std::cout<<"TODO: implement clipping. (The current code passes the triangle through without clipping them.)"<<std::endl;
    clip_triangle(state,in,face+1);
}

// Rasterize the triangle defined by the three vertices in the "in" array.  This
// function is responsible for rasterization, interpolation of data to
// fragments, calling the fragment shader, and z-buffering.
void rasterize_triangle(driver_state& state, const data_geometry* in[3])
{
    //std::cout<<"TODO: implement rasterization"<<std::endl;
    data_geometry out[3];
    data_vertex verts;

    for(unsigned int i = 0; i < 3; ++i){
        verts.data = in[i]->data;
        state.vertex_shader( verts, out[i], state.uniform_data);
    }

    float alpha, beta, gamma, inv_tArea;

    size_t width = state.image_width;
    size_t height = state.image_height;
    vec2 pixel_pos;
///*   
    vec2 a = vec2( (out[0].gl_Position[0] /out[0].gl_Position[3] + 1) * 0.5 * width ,
                   (out[0].gl_Position[1] /out[0].gl_Position[3] + 1) * 0.5 * height);

    vec2 b = vec2( (out[1].gl_Position[0] /out[1].gl_Position[3] + 1) * 0.5 * width ,
                   (out[1].gl_Position[1] /out[1].gl_Position[3] + 1) * 0.5 * height);

    vec2 c = vec2( (out[2].gl_Position[0] /out[2].gl_Position[3] + 1) * 0.5 * width ,
                   (out[2].gl_Position[1] /out[2].gl_Position[3] + 1) * 0.5 * height);
//*/
    inv_tArea = 1/t_Area(a,b,c);
    int walk = 0;
    for(size_t i = 0; i < height; ++i){
        for(size_t j = 0; j < width; ++j){
            
            pixel_pos = {float(j + 0.5), float(i + 0.5)};
            alpha = t_Area(pixel_pos, b, c) * inv_tArea;
            beta  = t_Area(a, pixel_pos, c) * inv_tArea;
            gamma = t_Area(a, b, pixel_pos) * inv_tArea;

            if(alpha >= 0 && beta >= 0 && gamma >= 0){
            //std::cout << "A B G = " << alpha << " " << beta << " " << gamma << std::endl;
                data_fragment dFrag;
                data_output dOutput;

                state.fragment_shader(dFrag, dOutput, state.uniform_data);
                //state.image_color[(j * width) + i] = make_pixel(255,255,255);
                state.image_color[walk] = make_pixel(255,255,255);
                //state.image_color[(i * width) + j] =make_pixel(dOutput.output_color[0], dOutput.output_color[1],dOutput.output_color[2]);
            }
            ++walk;
        }
    }
}

float t_Area(vec2 &a, vec2 &b, vec2 &c){
 return 0.5 * ((b[0]*c[1]-c[0]*b[1]) 
              -(a[0]*c[1]-c[0]*a[1])
              +(a[0]*b[1]-b[0]*a[1]));
}

