#include <iostream>
#include <GL/glew.h> 
#include <GLFW/glfw3.h>
#define CGAL_EIGEN3_ENABLED
#define GLM_FORCE_RADIANS
#include <glm/glm.hpp>
#include <glm/ext.hpp>
#include <iostream>
#include <string>
#include <assimp/Importer.hpp>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh_default_triangulation_3.h>
#include <CGAL/Complex_2_in_triangulation_3.h>
#include <CGAL/make_surface_mesh.h>
#include <CGAL/Implicit_surface_3.h>
#include <CGAL/IO/output_surface_facets_to_triangle_soup.h>
#include <CGAL/IO/output_surface_facets_to_polyhedron.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Extended_homogeneous.h>
#include <CGAL/Nef_polyhedron_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/IO/Nef_polyhedron_iostream_3.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/Polyhedron_items_with_id_3.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Boolean_set_operations_2.h>
#include <list>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Delaunay_mesher_2.h>
#include <CGAL/Delaunay_mesh_face_base_2.h>
#include <CGAL/Delaunay_mesh_size_criteria_2.h>
#include <CGAL/Triangulation_2.h>

// HalfedgeGraph adaptors for Polyhedron_3
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/boost/graph/properties_Polyhedron_3.h>
#include <CGAL/Null_matrix.h>
#include <CGAL/Surface_mesh_deformation.h>
#include <CGAL/Exact_integer.h>
#include <CGAL/Extended_homogeneous.h>
#include <CGAL/Nef_polyhedron_3.h>
#include <fstream>
#include <map>
#include <queue>

using namespace std;
using namespace Assimp;

const char* v_shader =
"#version 400\n"
"in vec3 vp;\n"
"uniform mat4 p;\n"
"uniform mat4 v;\n"
"uniform mat4 m;\n"
"out vec3 vp_color;\n"
"void main() {\n"
"   gl_Position = p * v * m * vec4(vp, 1.0);\n"
"   vp_color = vp;\n"
"}\n";

const char* f_shader =
"#version 400\n"
"out vec4 frag_color;\n"
"in vec3 vp_color;\n"
"void main() {\n"
"   frag_color = vec4(vp_color, 0.5);\n"
"}\n";

GLfloat points[] = {
	0.0f, 0.5f, 0.0f, //0

	-0.5f, -0.5f, -0.5f,  //1
	0.5f, -0.5f, -0.5f,  //2
	0.5f, -0.5f,  0.5f,  //3
	-0.5f, -0.5f,  0.5f,  //4
};

GLushort indexs[] = {
	0, 1, 2, //ÕýÃæ
	0, 2, 3, //ÓÒ±ß
	0, 3, 4, //±³Ãæ
	0, 4, 1, //×ó±ß
	1, 4, 3, 1, 3, 2 //ÏÂÃæ
};

GLuint vbo = 0;
GLuint vao = 0;
GLuint ibo = 0;
GLuint program = 0;

void error_out(int err, const char* desc) {
	cout << "Error: " << err << " || " << desc << endl;
};

// default triangulation for Surface_mesher
typedef CGAL::Surface_mesh_default_triangulation_3 Tr;
// c2t3
typedef CGAL::Complex_2_in_triangulation_3<Tr> C2t3;
typedef Tr::Geom_traits GT;
typedef GT::Sphere_3 Sphere_3;
typedef GT::Point_3 Point_3;
typedef GT::FT FT;
typedef FT (*Function)(Point_3);
typedef CGAL::Implicit_surface_3<GT, Function> Surface_3;

//typedef CGAL::Simple_cartesian<CGAL::Gmpz> Kernel;
typedef CGAL::Epick Kernel;
//typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel2;
typedef CGAL::Homogeneous<CGAL::Epick>  Kernel2;
typedef CGAL::Polyhedron_3<Kernel,CGAL::Polyhedron_items_with_id_3>  Polyhedron;
typedef CGAL::Polyhedron_3<Kernel2>  Polyhedron2;
typedef boost::graph_traits<Polyhedron>::vertex_descriptor    vertex_descriptor;
typedef boost::graph_traits<Polyhedron>::vertex_iterator        vertex_iterator;
typedef boost::graph_traits<Polyhedron>::halfedge_descriptor halfedge_descriptor;
typedef boost::graph_traits<Polyhedron>::out_edge_iterator    out_edge_iterator;
typedef CGAL::Nef_polyhedron_3<Kernel2>  Nef_polyhedron;
typedef CGAL::Cartesian_converter<Kernel2, Kernel> converter;

FT sphere_function (Point_3 p) {
    const FT x2=p.x()*p.x(), y2=p.y()*p.y(), z2=p.z()*p.z();
    return x2+y2+z2-1;
}

// Can be used to convert polyhedron from exact to inexact and vice-versa
template <class Polyhedron_input,
class Polyhedron_output>
struct Copy_polyhedron_to
: public CGAL::Modifier_base<typename Polyhedron_output::HalfedgeDS>
{
    Copy_polyhedron_to(const Polyhedron_input& in_poly)
    : in_poly(in_poly) {}
    
    void operator()(typename Polyhedron_output::HalfedgeDS& out_hds)
    {
        typedef typename Polyhedron_output::HalfedgeDS Output_HDS;
        typedef typename Polyhedron_input::HalfedgeDS Input_HDS;
        
        CGAL::Polyhedron_incremental_builder_3<Output_HDS> builder(out_hds);
        
        typedef typename Polyhedron_input::Vertex_const_iterator Vertex_const_iterator;
        typedef typename Polyhedron_input::Facet_const_iterator  Facet_const_iterator;
        typedef typename Polyhedron_input::Halfedge_around_facet_const_circulator HFCC;
        
        builder.begin_surface(in_poly.size_of_vertices(),
                              in_poly.size_of_facets(),
                              in_poly.size_of_halfedges());
        
        for(Vertex_const_iterator
            vi = in_poly.vertices_begin(), end = in_poly.vertices_end();
            vi != end ; ++vi)
        {
            typename Polyhedron_output::Point_3 p(::CGAL::to_double( vi->point().x()),
                                                  ::CGAL::to_double( vi->point().y()),
                                                  ::CGAL::to_double( vi->point().z()));
            builder.add_vertex(p);
        }
        
        typedef CGAL::Inverse_index<Vertex_const_iterator> Index;
        Index index( in_poly.vertices_begin(), in_poly.vertices_end());
        
        for(Facet_const_iterator
            fi = in_poly.facets_begin(), end = in_poly.facets_end();
            fi != end; ++fi)
        {
            HFCC hc = fi->facet_begin();
            HFCC hc_end = hc;
            builder.begin_facet ();
            do {
                builder.add_vertex_to_facet(index[hc->vertex()]);
                ++hc;
            } while( hc != hc_end);
            builder.end_facet();
        }
        builder.end_surface();
    } // end operator()(..)
private:
    const Polyhedron_input& in_poly;
}; // end Copy_polyhedron_to<>

template <class Poly_B, class Poly_A>
void poly_copy(Poly_B& poly_b, const Poly_A& poly_a)
{
    poly_b.clear();
    Copy_polyhedron_to<Poly_A, Poly_B> modifier(poly_a);
    poly_b.delegate(modifier);
}

// Collect the vertices which are at distance less or equal to k
// from the vertex v in the graph of vertices connected by the edges of P
std::vector<vertex_descriptor> extract_k_ring(const Polyhedron &P, vertex_descriptor v, int k)
{
    std::map<vertex_descriptor, int>  D;
    std::vector<vertex_descriptor>    Q;
    Q.push_back(v); D[v] = 0;
    std::size_t current_index = 0;
    int dist_v;
    while( current_index < Q.size() && (dist_v = D[ Q[current_index] ]) < k ) {
        v = Q[current_index++];
        out_edge_iterator e, e_end;
        for(boost::tie(e, e_end) = out_edges(v, P); e != e_end; e++)
        {
            halfedge_descriptor he = halfedge(*e, P);
            vertex_descriptor new_v = target(he, P);
            if(D.insert(std::make_pair(new_v, dist_v + 1)).second) {
                Q.push_back(new_v);
            }
        }
    }
    return Q;
}


vector<Point_3> get_tr() {
	Polyhedron poly1, poly2;
    Tr tr;            // 3D-Delaunay triangulation
    C2t3 c2t3 (tr);   // 2D-complex in 3D-Delaunay triangulation
    // defining the surface
    std::cout << "1 number of points: " << tr.number_of_vertices() << "\n";
    Surface_3 surface(sphere_function,             // pointer to function
                      Sphere_3(CGAL::ORIGIN, 2.)); // bounding sphere
    // Note that "2." above is the *squared* radius of the bounding sphere!
    // defining meshing criteria
    std::cout << "2 number of points: " << tr.number_of_vertices() << "\n";
    CGAL::Surface_mesh_default_criteria_3<Tr> criteria(30.,  // angular bound
                                                       0.1,  // radius bound
                                                       0.1); // distance bound
    // meshing surface
    std::cout << "3 number of points: " << tr.number_of_vertices() << "\n";
    CGAL::make_surface_mesh(c2t3, surface, criteria, CGAL::Non_manifold_tag());
    std::cout << "Final number of points: " << tr.number_of_vertices() << "\n";
    
    Polyhedron poly;
    CGAL::output_surface_facets_to_polyhedron(c2t3, poly);
    CGAL::Surface_mesh_deformation<Polyhedron> deform(poly);

	// KENT WARN : 这一步很重要，没有会挂掉
	set_halfedgeds_items_id(poly);

	// Select and insert the vertices of the region of interest
	vertex_iterator vb, ve;
    boost::tie(vb,ve) = CGAL::vertices(poly);
	std::vector<vertex_descriptor> roi = extract_k_ring(poly, *CGAL::cpp11::next(vb, 0), poly.size_of_vertices());
	deform.insert_roi_vertices(roi.begin(), roi.end());

	// Select and insert the control vertices
	std::vector<vertex_descriptor> cvertices_1 = extract_k_ring(poly, *CGAL::cpp11::next(vb, 0), 20);
	deform.insert_control_vertices(cvertices_1.begin(), cvertices_1.end());

    deform.translate(cvertices_1.begin(), cvertices_1.end(), Eigen::Vector3d(0,0.4,0));
	deform.deform();

	// 构造上面那个球形完成
	poly1 = poly;

	// 构造下面那个球形开始
	CGAL::output_surface_facets_to_polyhedron(c2t3, poly2);
	CGAL::Surface_mesh_deformation<Polyhedron> deform2(poly2);

	// KENT WARN : 这一步很重要，没有会挂掉
	set_halfedgeds_items_id(poly2);

	// Select and insert the vertices of the region of interest
    boost::tie(vb,ve) = CGAL::vertices(poly2);
	roi = extract_k_ring(poly2, *CGAL::cpp11::next(vb, 0), poly2.size_of_vertices());
	deform2.insert_roi_vertices(roi.begin(), roi.end());

	// Select and insert the control vertices
	cvertices_1 = extract_k_ring(poly2, *CGAL::cpp11::next(vb, 0), 20);
	deform2.insert_control_vertices(cvertices_1.begin(), cvertices_1.end());

    deform2.translate(cvertices_1.begin(), cvertices_1.end(), Eigen::Vector3d(0,-0.4,0));
	deform2.deform();
	// 构造下面那个球形完成

	// 进行布尔操作
    Polyhedron2 p1;
    poly_copy(p1, poly1);
	Nef_polyhedron ball1(p1);

	ball1.convert_to_polyhedron(p1);

	// OUTPUT points
    vector<Point_3> pts;
	for (auto it = poly.facets_begin(); it != poly.facets_end(); it++) {
		auto it2 = it->facet_begin();
		auto it3 = it2;
		CGAL_For_all(it2, it3) {
			pts.push_back(Point_3(it2->vertex()->point().x(),
								  it2->vertex()->point().y(),
								  it2->vertex()->point().z()));
		}
	}

    cout << pts.size() << endl;

    return pts;
}


int main(int argc, const char * argv[]) {
	int w = 1024;
	int h = 768;
    
	if (!glfwInit()) {
		cout << "error to init glfw" << endl;
	}
    
	glfwWindowHint(GLFW_SAMPLES, 4);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
	glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
	glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);

	glfwSetErrorCallback(error_out);

	GLFWwindow *win = glfwCreateWindow(w, h, "this is glfw window", nullptr, nullptr);
	if (!win) {
		cout << "error to create window" << endl;
		glfwTerminate();
		return -1;
	}

	glfwMakeContextCurrent(win);
	glfwSwapInterval(1);

	glewExperimental = GL_TRUE;
	GLenum err = glewInit();
	if (GLEW_OK != err)
	{

		fprintf(stderr, "Error: %s/n", glewGetErrorString(err));
	}
	fprintf(stdout, "Status: Using GLEW %s\n", glewGetString(GLEW_VERSION));

	const GLubyte* name = glGetString(GL_VENDOR); //·µ»Ø¸ºÔðµ±Ç°OpenGLÊµÏÖ³§ÉÌµÄÃû×Ö
	const GLubyte* biaoshifu = glGetString(GL_RENDERER); //·µ»ØÒ»¸öäÖÈ¾Æ÷±êÊ¶·û£¬Í¨³£ÊÇ¸öÓ²¼þÆ½Ì¨
	const GLubyte* OpenGLVersion = glGetString(GL_VERSION); //·µ»Øµ±Ç°OpenGLÊµÏÖµÄ°æ±¾ºÅ
	printf("OpenGL%s\n", name);
	printf("%s\n", biaoshifu);
	printf("OOpenGL%s\n", OpenGLVersion);

	//create shaders
	/*GLuint vs = glCreateShader();*/
	GLuint vs = glCreateShader(GL_VERTEX_SHADER);

	glShaderSource(vs, 1, &v_shader, nullptr);
	glCompileShader(vs);
	GLuint fs = glCreateShader(GL_FRAGMENT_SHADER);
	glShaderSource(fs, 1, &f_shader, nullptr);
	glCompileShader(fs);

	//shader error
	    GLint *status = new GLint;
	    glGetShaderiv(fs, GL_COMPILE_STATUS, status);
	    if(*status == GL_FALSE){
	        GLint length;
	        glGetShaderiv(fs, GL_INFO_LOG_LENGTH, &length);
	        string message("Failed to create ");
	        message += (GL_FRAGMENT_SHADER == GL_VERTEX_SHADER ? "vertex shade:\n" : "fragment shader:\n");
	        
	        char *detail = new char[length];
	        glGetShaderInfoLog(fs, length, &length, detail);
	        message += detail;
	        delete[] detail;
	        
	        glDeleteShader(fs);
	        
	        cout<<message<<endl;
	    }

	//create program
	program = glCreateProgram();
	glAttachShader(program, vs);
	glAttachShader(program, fs);
	glLinkProgram(program);

    vector<Point_3> pts = get_tr();
    int pts_size = pts.size() * 3;
    GLfloat  * pts_data = (GLfloat *)malloc(sizeof(GLfloat) * pts_size);
    int i = 0;
    for (auto it = pts.begin(); it != pts.end(); it++) {
        pts_data[3 * i + 0] = it->x();
        pts_data[3 * i + 1] = it->y();
        pts_data[3 * i + 2] = it->z();
        cout << *it << endl;
        i++;
    }
    
    //set vbo
	glGenBuffers(1, &vbo);
	glBindBuffer(GL_ARRAY_BUFFER, vbo);
	glBufferData(GL_ARRAY_BUFFER, (sizeof(GLfloat) * pts_size), pts_data, GL_STATIC_DRAW);
    
	//set ibo
	glGenBuffers(1, &ibo);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ibo);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(indexs), indexs, GL_STATIC_DRAW);

	//use vbo
	glGenVertexArrays(1, &vao);
	glBindVertexArray(vao);
	glEnableVertexAttribArray(0);
	glBindBuffer(GL_ARRAY_BUFFER, vbo);
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, nullptr);

	GLfloat angle = 0.0f;

	while (!glfwWindowShouldClose(win)) {
		angle += 0.05f;

		glfwGetFramebufferSize(win, &w, &h);
		glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
		glEnable(GL_DEPTH_TEST);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		//ÊÓ½Ç¿í£¨½Ç¶È£©       ¿í¸ß±È    ½üÇÐÃæ  Ô¶ÇÐÃæ
		glm::mat4 proj = glm::perspective(glm::radians(60.0f), (float)w / h, 0.3f, 500.0f);
		//ÑÛ¾¦Î»ÖÃ      ¿´µÄµãµÄÎ»ÖÃ       Í·¶¥µÄÏòÁ¿
		glm::mat4 camera = glm::lookAt(glm::vec3(0, 0, 3), glm::vec3(0, 0, 0), glm::vec3(0, 1, 0));
		glUniformMatrix4fv(glGetUniformLocation(program, "p"), 1, GL_FALSE, glm::value_ptr(proj));
		glUniformMatrix4fv(glGetUniformLocation(program, "v"), 1, GL_FALSE, glm::value_ptr(camera));
		glUniformMatrix4fv(glGetUniformLocation(program, "m"), 1, GL_FALSE, glm::value_ptr(glm::rotate(glm::radians(angle), glm::vec3(0.0f, 1.0f, 0.0f))));

		glUseProgram(program);
		glBindVertexArray(vao);

		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ibo);
        glDrawArrays(GL_TRIANGLES, 0, pts_size);
//		glDrawElements(GL_TRIANGLES, 18, GL_UNSIGNED_SHORT, nullptr);
		glfwPollEvents();
		glfwSwapBuffers(win);
	}

	glfwDestroyWindow(win);

	new Importer;
};


