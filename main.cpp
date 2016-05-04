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
// HalfedgeGraph adaptors for Polyhedron_3
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/boost/graph/properties_Polyhedron_3.h>
#include <CGAL/Null_matrix.h>
#include <CGAL/Surface_mesh_deformation.h>
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

typedef CGAL::Epick Kernel;
typedef CGAL::Polyhedron_3<Kernel,CGAL::Polyhedron_items_with_id_3>  Polyhedron;
typedef boost::graph_traits<Polyhedron>::vertex_descriptor    vertex_descriptor;
typedef boost::graph_traits<Polyhedron>::vertex_iterator        vertex_iterator;
typedef boost::graph_traits<Polyhedron>::halfedge_descriptor halfedge_descriptor;
typedef boost::graph_traits<Polyhedron>::out_edge_iterator    out_edge_iterator;

FT sphere_function (Point_3 p) {
    const FT x2=p.x()*p.x(), y2=p.y()*p.y(), z2=p.z()*p.z();
    return x2+y2+z2-1;
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
    
    vertex_iterator vb, ve;
    boost::tie(vb,ve) = CGAL::vertices(poly);

    std::vector<vertex_descriptor> cvertices_1 = extract_k_ring(poly, *CGAL::cpp11::next(vb, 0), 20);
    std::vector<vertex_descriptor> cvertices_2 = extract_k_ring(poly, *CGAL::cpp11::next(vb, 97), 20);
    
//    deform.insert_control_vertices(cvertices_1.begin(), cvertices_1.end());
//
//    deform.translate(cvertices_1.begin(), cvertices_1.end(), Eigen::Vector3d(0,0.3,0));
    
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


