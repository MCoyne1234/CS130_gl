#include "driver_state.h"
#include <cstring>
#include <vector>

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
    state.image_depth = new float[area];//0;
    for (size_t i = 0; i < area;++i) {
        state.image_depth[i] = -200.0;
    }

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

///*
        //int fpv = state.floats_per_vertex;
        data_geometry dg_arr[3];
        const data_geometry * p_dg_arr[3] = {&dg_arr[0], &dg_arr[1], &dg_arr[2]};
        for(int vert = 0 ;vert < state.num_vertices * state.floats_per_vertex ; vert += 3*state.floats_per_vertex){
            for(int v = 0; v < 3; ++v){
                dg_arr[v].data = state.vertex_data+vert+v*state.floats_per_vertex;
            }
            clip_triangle(state, p_dg_arr, 0);
            //rasterize_triangle(state, p_dg_arr);
        }
//*/

        /*
            data_geometry dg_arr[3];
            const data_geometry * p_dg_arr[3] = {&dg_arr[0], &dg_arr[1], &dg_arr[2]};
            int incr = 3 * state.floats_per_vertex;
            for(int i = 0; i < state.num_vertices*state.floats_per_vertex; i += incr){
                std::cout << "YO " << state.num_vertices << std::endl;
                    dg_arr[0].data = (state.vertex_data+i);
                    dg_arr[1].data = (state.vertex_data+i+state.floats_per_vertex);
                    dg_arr[2].data = (state.vertex_data+i+2*state.floats_per_vertex);
                    rasterize_triangle(state, p_dg_arr);
            }
         */
            }
            break;
        case render_type::indexed :
            {
            //std::cout << "indexed" << std::endl;
            data_geometry dg_index[3];
            const data_geometry * p_dg_index[3] = {&dg_index[0],&dg_index[1],&dg_index[2]};
            size_t max =  3 * state.num_triangles;

            for(size_t i = 0; i < max; i += 3){
                dg_index[0].data = state.vertex_data + *(state.index_data +i) * 3;//state.floats_per_vertex;
                dg_index[1].data = state.vertex_data + *(state.index_data +i+1) * 3;//state.floats_per_vertex;
                dg_index[2].data = state.vertex_data + *(state.index_data +i+2) * 3;//state.floats_per_vertex;
                clip_triangle(state, p_dg_index, 0);
                //rasterize_triangle(state, p_dg_index);
            }
            break;
            }
        case render_type::fan :
            {
            //std::cout << "fan" << std::endl;
            data_geometry dg_center[3];
            const data_geometry * p_dg_center[3] = {&dg_center[0],&dg_center[1],&dg_center[2]};
            size_t max = state.num_vertices * state.floats_per_vertex;
            size_t incr = 3 * state.floats_per_vertex;

            dg_center[0].data = state.vertex_data;
            dg_center[1].data = state.vertex_data + state.floats_per_vertex;
            dg_center[2].data = state.vertex_data + 2*state.floats_per_vertex;
            clip_triangle(state, p_dg_center, 0);
            //rasterize_triangle(state, p_dg_center);
            for(size_t i = incr; i < max; i += state.floats_per_vertex){
                dg_center[1].data = dg_center[2].data;
                dg_center[2].data = state.vertex_data + i;

                clip_triangle(state, p_dg_center, 0);
                //rasterize_triangle(state, p_dg_center);
            }

            break;
            }
        case render_type::strip :
            {
            //std::cout << "strip" << std::endl;
            data_geometry dg_arr[3];
            const data_geometry * p_dg_arr[3] = {&dg_arr[0],&dg_arr[1],&dg_arr[2]};
            size_t max = state.num_vertices * state.floats_per_vertex;
            size_t incr = 3 * state.floats_per_vertex;

            dg_arr[0].data = state.vertex_data;
            dg_arr[1].data = state.vertex_data + state.floats_per_vertex;
            dg_arr[2].data = state.vertex_data + 2*state.floats_per_vertex;
            clip_triangle(state, p_dg_arr, 0);

            for(size_t i = incr; i < max; i += state.floats_per_vertex){
                dg_arr[0].data = dg_arr[1].data;
                dg_arr[1].data = dg_arr[2].data;
                dg_arr[2].data = state.vertex_data + i;
                clip_triangle(state, p_dg_arr, 0);
            }
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
void clip_triangle(driver_state& state, const data_geometry* in[3], int face)
{
    data_vertex verts;
    data_geometry out[3];
    const data_geometry * p_out[3] = {&out[0] ,&out[1], &out[2]};

        for(unsigned int i = 0; i < 3; ++i){
            verts.data = in[i]->data;
            out[i].data = verts.data;
            out[i].gl_Position = in[i]->gl_Position;
    if (face == 0){
            state.vertex_shader(verts, out[i], state.uniform_data);
        }
    }
    //rasterize_triangle(state, p_out);
    //return;
    /*
    vec3 a = {out[0].data[0] , out[0].data[1] , out[0].data[2]};
    vec3 b = {out[1].data[0] , out[1].data[1] , out[1].data[2]};
    vec3 c = {out[2].data[0] , out[2].data[1] , out[2].data[2]};
    */
    ///*
    vec3 a = {out[0].gl_Position[0] , out[0].gl_Position[1] , out[0].gl_Position[2]};
    vec3 b = {out[1].gl_Position[0] , out[1].gl_Position[1] , out[1].gl_Position[2]};
    vec3 c = {out[2].gl_Position[0] , out[2].gl_Position[1] , out[2].gl_Position[2]};
    //*/
    vec3 plane_point;
    vec3 plane_norm;

    std::vector<vec3> clipped_in, clipped_out;
    vec3 clipped_v1, clipped_v2, u,v;//u1, v1, u2, v2;

    if(face==6){
       //std::cout << "raster.raster.raster" <<  std::endl;
        rasterize_triangle(state, p_out);
        return;
    }else if(face==5){
        plane_point ={0,1,0};
        plane_norm = {0,-1,0};

    }else if(face==4){
        plane_point ={0,-1,0};
        plane_norm = {0,1,0};

    }else if(face==3){
        plane_point ={1,0,0};
        plane_norm = {-1,0,0};

    }else if(face==2){
        plane_point ={-1,0,0};
        plane_norm = {1,0,0};

    }else if(face==1){
        plane_point ={0,0,-1};
        plane_norm = {0,0,1};

    }else if(face==0){
       plane_point ={0,0,1};
       plane_norm = {0,0,-1};

    }
       /*
       float check_a = dot(plane_norm,(plane_point - a));
       float check_b = dot(plane_norm,(plane_point - b));
       float check_c = dot(plane_norm,(plane_point - c));

       std::cout << "check a " << check_a << std::endl;
       std::cout << "check b " << check_b << std::endl;
       std::cout << "check c " << check_c << std::endl;
       std::cout << "W = " << out[vert_index].gl_Position[3] << std::endl << std::endl;
       */

        std::vector<data_geometry> holder;
        std::vector<data_geometry> goner;
        if(dot(plane_norm,(plane_point - a)) <= out[0].gl_Position[3]){
            clipped_in.push_back(a);
            holder.push_back(out[0]);
        }else {
            clipped_out.push_back(a);
            goner.push_back(out[0]);
        }

        if(dot(plane_norm,(plane_point - b)) <= out[1].gl_Position[3]){
            clipped_in.push_back(b);
            holder.push_back(out[1]);
        }else {
            clipped_out.push_back(b);
            goner.push_back(out[1]);
        }

        if(dot(plane_norm,(plane_point - c)) <= out[2].gl_Position[3]){
            clipped_in.push_back(c);
            holder.push_back(out[2]);
        }else {
            clipped_out.push_back(c);
            goner.push_back(out[2]);
        }

        float D,N;
        data_geometry new_out_1[3];
        const data_geometry * p_new_out_1[3] = {&new_out_1[0] ,&new_out_1[1] ,&new_out_1[2]} ;
        data_geometry new_out_2[3];
        const data_geometry * p_new_out_2[3] = {&new_out_2[0] ,&new_out_2[1] ,&new_out_2[2]} ;

        //float new_data_A[state.floats_per_vertex] = {0,0,0};
        //float new_data_B[state.floats_per_vertex] = {0,0,0};

        //float * p_A= &new_data_A[0];
        //float * p_B= &new_data_B[0];

        if(clipped_in.size() == 0){ return;}// std::cout << "DERP" << std::endl; }
        if(clipped_in.size() == 1){

            //std::cout << "TWO POINTS OUT " << std::endl;

            new_out_1[0].data = holder[0].data;
            new_out_1[0].gl_Position = holder[0].gl_Position;

            // GET FIRST NEW POINT
            u = clipped_out[0] - clipped_in[0];
            v = clipped_in[0] - plane_point;
            D = dot(plane_norm,u);
            N = dot(plane_norm, v);

            float sI = N / D;
            clipped_v1 = clipped_in[0] + sI * u;

            new_out_1[1].data = goner[0].data;
            new_out_1[1].gl_Position = goner[0].gl_Position;

            new_out_1[1].gl_Position[0] = clipped_v1[0];
            new_out_1[1].gl_Position[1] = clipped_v1[1];
            new_out_1[1].gl_Position[2] = clipped_v1[2];


            // GET SECOND NEW POINT
            u = clipped_out[1] - clipped_in[0];
            v = clipped_in[0] - plane_point;
            D = dot(plane_norm, u);
            N = dot(plane_norm, v);
            sI = N / D;

            clipped_v2 = clipped_in[0] + sI * u;

            new_out_1[2].data = goner[1].data;
            new_out_1[2].gl_Position = goner[1].gl_Position;

            new_out_1[2].gl_Position[0] = clipped_v2[0];
            new_out_1[2].gl_Position[1] = clipped_v2[1];
            new_out_1[2].gl_Position[2] = clipped_v2[2];

            clip_triangle(state, p_new_out_1, face+1);
        }
        if(clipped_in.size() == 2){
        //std::cout << "ONE POINT OUT " << std::endl;
/*
            u = clipped_out[0] - clipped_in[0];
            v = clipped_in[0] - plane_point;
            D = dot(plane_norm,u);
            N = dot(plane_norm, v);
            float sI = N / D;

            clipped_v1 = clipped_in[0] + sI * u;

            u = clipped_out[0] - clipped_in[1];
            v = clipped_in[1] - plane_point;
            D = dot(plane_norm,u);
            N = dot(plane_norm, v);
            sI = N / D;

            clipped_v2 = clipped_in[1] + sI * u;

            for(int i = 0; i < state.floats_per_vertex; ++i){
                //std::cout << "DERKA-DER " << face<< std::endl;
                //std::cout << "DER-DER " << (goner[0]).data[0] << std::endl;
                new_data_A[i] = goner[0].data[i];
                new_data_B[i] = goner[0].data[i];
            }

            new_data_A[0] *= clipped_v1[0];
            new_data_A[1] *= clipped_v1[1];
            new_data_A[2] *= clipped_v1[2];

            new_data_B[0] *=  clipped_v2[0];
            new_data_B[1] *=  clipped_v2[1];
            new_data_B[2] *=  clipped_v2[2];

            // FIRST NEW TRiANGLE
            new_out_1[0].data = holder[0].data;
            new_out_1[0].gl_Position = holder[0].gl_Position;

            new_out_1[1].data = p_A;
            //new_out_1[1].data = goner[0].data;

            new_out_1[1].gl_Position = goner[0].gl_Position;
            new_out_1[1].gl_Position[0] = clipped_v1[0];
            new_out_1[1].gl_Position[1] = clipped_v1[1];
            new_out_1[1].gl_Position[2] = clipped_v1[2];

            new_out_1[2].data = p_B;
           // new_out_1[2].data = goner[0].data;
            new_out_1[2].gl_Position = goner[0].gl_Position;

            new_out_1[2].gl_Position[0] = clipped_v2[0];
            new_out_1[2].gl_Position[1] = clipped_v2[1];
            new_out_1[2].gl_Position[2] = clipped_v2[2];


            // SECOND NEW TRIANGLE
            new_out_2[0].data = holder[0].data;
            new_out_2[0].gl_Position = holder[0].gl_Position;

            new_out_2[1].data = holder[1].data;
            new_out_2[1].gl_Position = holder[1].gl_Position;

            new_out_1[2].data = p_B;
            //new_out_2[2].data = goner[0].data;

            new_out_2[2].gl_Position = goner[0].gl_Position;
            new_out_2[2].gl_Position[0] = clipped_v2[0];
            new_out_2[2].gl_Position[1] = clipped_v2[1];
            new_out_2[2].gl_Position[2] = clipped_v2[2];

            clip_triangle(state, p_new_out_1, face+1);
            clip_triangle(state, p_new_out_2, face+1);
            return;
       */
        }
        //dot(plane_norm,(plane_point - a));

       clip_triangle(state, p_out, face+1);

    //std::cout<<"TODO: implement clipping. (The current code passes the triangle through without clipping them.)"<<std::endl;
    //clip_triangle(state,in,face+1);
}

// Rasterize the triangle defined by the three vertices in the "in" array.  This
// function is responsible for rasterization, interpolation of data to
// fragments, calling the fragment shader, and z-buffering.
void rasterize_triangle(driver_state& state, const data_geometry* in[3])
{
/*
    //std::cout<<"TODO: implement rasterization"<<std::endl;
    std::cout<<"a: "<< in[0]->gl_Position<<std::endl;
    std::cout<<"b: "<< in[1]->gl_Position<<std::endl;
    std::cout<<"c: "<< in[2]->gl_Position<<std::endl;
*/

    float alpha, beta, gamma, inv_tArea;

    size_t width = state.image_width;
    size_t height = state.image_height;
 ///*
    vec2 pixel_pos;
    vec3 pixel_pos_3d;

    vec2 a = vec2( (in[0]->gl_Position[0] /in[0]->gl_Position[3] + 1) * 0.5 * width ,
                   (in[0]->gl_Position[1] /in[0]->gl_Position[3] + 1) * 0.5 * height);

    vec2 b = vec2( (in[1]->gl_Position[0] /in[1]->gl_Position[3] + 1) * 0.5 * width ,
                   (in[1]->gl_Position[1] /in[1]->gl_Position[3] + 1) * 0.5 * height);

    vec2 c = vec2( (in[2]->gl_Position[0] /in[2]->gl_Position[3] + 1) * 0.5 * width ,
                   (in[2]->gl_Position[1] /in[2]->gl_Position[3] + 1) * 0.5 * height);
//*/

    inv_tArea = 1/t_Area(a,b,c);
    int walk = 0;
    for(size_t i = 0; i < height; ++i){
        for(size_t j = 0; j < width; ++j){
            
            pixel_pos = {float(j + 0.5), float(i + 0.5)};
            alpha = t_Area(pixel_pos, b, c) * inv_tArea;
            beta  = t_Area(a, pixel_pos, c) * inv_tArea;
            gamma = t_Area(a, b, pixel_pos) * inv_tArea;

            pixel_pos_3d[0] = pixel_pos[0];
            pixel_pos_3d[1] = pixel_pos[1];

            //pixel_pos_3d[2] = (in[0]->data[2]*alpha + in[1]->data[2]*beta + in[2]->data[2]*gamma);
            //std::cout<<"data: "<<  pixel_pos_3d[2]<<std::endl;

            //pixel_pos_3d[2] = -(in[0]->gl_Position[2]*alpha + in[1]->gl_Position[2]*beta + in[2]->gl_Position[2]*gamma);
            //std::cout<<"gl_pos: "<< pixel_pos_3d[2]<<std::endl;
            pixel_pos_3d[2] = -(in[0]->gl_Position[2]/in[0]->gl_Position[3]*alpha + in[1]->gl_Position[2]/in[1]->gl_Position[3]*beta + in[2]->gl_Position[2]/in[2]->gl_Position[3]*gamma);

            if(alpha >= 0 && beta >= 0 && gamma >= 0){
            //std::cout << "A B G = " << alpha << " " << beta << " " << gamma << std::endl;
                float d_out[state.floats_per_vertex];
                //data_fragment dFrag;
                data_output dOutput;


                for(int k = 0; k < state.floats_per_vertex; ++k){
                    if(state.interp_rules[k] == interp_type::flat){
                        // then data_out[i] should be the float from data from the first vertex.
                        d_out[k]= in[0]->data[k];
                    }
                    else if(state.interp_rules[k] == interp_type::noperspective){
                        // then data_out[i] should be the interpolation of data_a[i], data_b[i] and data_c[i] using the barycentric coordinates.
                        d_out[k]= (alpha * in[0]->data[k]) + (beta * in[1]->data[k]) + (gamma * in[2]->data[k]);
                    }
                    else if(state.interp_rules[k] == interp_type::smooth){
                        // then data_out[i] should be the interpolation of data_a[i], data_b[i] and data_c[i]
                        // using the barycentric coordinates with perspective correction
                        //d_out[k]= (alpha * in[0]->data[k]/in[0]->gl_Position[3]) + (beta * in[1]->data[k]/in[1]->gl_Position[3]) + (gamma * in[2]->data[k]/in[2]->gl_Position[3]);

                        d_out[k]= ((alpha * in[0]->data[k]/in[0]->gl_Position[3]) + (beta * in[1]->data[k]/in[1]->gl_Position[3]) + (gamma * in[2]->data[k]/in[2]->gl_Position[3]))
                                /((alpha/in[0]->gl_Position[3]) + (beta/in[1]->gl_Position[3]) + (gamma/in[2]->gl_Position[3]));

                        // super wrong, but cool;
                        /*d_out[k]= ((alpha * in[0]->data[k]) + (beta * in[1]->data[k]) + (gamma * in[2]->data[k]))
                        * ( in[0]->gl_Position[3] * in[1]->gl_Position[3] * in[2]->gl_Position[3]);*/
                    }
                }

              // float hmm = (in[0]->gl_Position[]

               //if(pixel_pos_3d[2] > state.image_depth[walk] ){
               if(pixel_pos_3d[2] > state.image_depth[walk] && pixel_pos_3d[2] >= -1 && pixel_pos_3d[2] <= 1  ){
                    state.image_depth[walk] = pixel_pos_3d[2];
                    state.fragment_shader({d_out}, dOutput, state.uniform_data);
                    state.image_color[walk] = make_pixel(255*dOutput.output_color[0],255*dOutput.output_color[1],255*dOutput.output_color[2]);
                }
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

