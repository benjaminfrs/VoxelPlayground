#include <iostream>
#include <cmath>

#define GL_SILENCE_DEPRECATION

// Without this gl.h gets included instead of gl3.h
#define GLFW_INCLUDE_NONE
#include <GLFW/glfw3.h>

// For includes related to OpenGL, make sure their are included after glfw3.h
#include <OpenGL/gl3.h>

#include "shaders_class.hpp"
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

void keyCallback(GLFWwindow* window, int key, int scancode, int action, int mods)
{
    if (key == GLFW_KEY_ESCAPE && action == GLFW_PRESS)
        glfwSetWindowShouldClose(window, GLFW_TRUE);
}

void frameBufferResizeCallback(GLFWwindow* window, int width, int height){
   glViewport(0, 0, width, height);
}

int main(void)
{
    //Init and configure
    glfwInit();
    GLFWwindow* window;
    // Without these two hints, nothing above OpenGL version 2.1 is supported
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GLFW_TRUE);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);

    // Create a windowed mode window and its OpenGL context
    window = glfwCreateWindow(640, 480, "Running OpenGL on Mac", NULL, NULL);
    if (!window)
    {
        glfwTerminate();
        return -1;
    }
    // Make the window's context current
    glfwMakeContextCurrent(window);
    // Set callback fro framebuffer
    glfwSetFramebufferSizeCallback(window, frameBufferResizeCallback);
    // Set callback for window
    glfwSetKeyCallback(window, keyCallback);
    //------------------------------------Init and Configure

    //Build and compile shaders
    //----------------------------------------------------------
    Shader ourShader("vertex_shader.vs", "fragment_shader.fs"); 

    //Set up vertex data, buffers, and vertex attributes
    // Vertex data and buffer
		float vertices[] = {
		     0.5f,  0.5f, 0.0f,  // top right
		     0.5f, -0.5f, 0.0f,  // bottom right
		    -0.5f, -0.5f, 0.0f,  // bottom left
		    -0.5f,  0.5f, 0.0f   // top left 
		};
		unsigned int indices[] = {  // note that we start from 0!
		    0, 1, 3,   // first triangle
		    1, 2, 3    // second triangle
		};  
    /* Vertex array object*/
    unsigned int VAO, VBO, EBO;
    glGenVertexArrays(1, &VAO); 
    glGenBuffers(1, &VBO); 
    glGenBuffers(1, &EBO); 

    // bind the Vertex Array Object first, then bind and set vertex buffer(s), and then configure vertex attributes(s).
    glBindVertexArray(VAO);

		glBindBuffer(GL_ARRAY_BUFFER, VBO);
		glBufferData(GL_ARRAY_BUFFER, sizeof(vertices), vertices, GL_STATIC_DRAW);


		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
		glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(indices), indices, GL_STATIC_DRAW);


		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
		glDrawElements(GL_TRIANGLES, 6, GL_UNSIGNED_INT, 0);

    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);  

    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glBindVertexArray(0);

    //testing wireframe
    //glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

    // Render loop 
		unsigned int transformLoc;
		glm::mat4 trans = glm::mat4(1.0f);
    while (!glfwWindowShouldClose(window))
    {
        // Resize the viewport
        int width, height;
        glfwGetFramebufferSize(window, &width, &height);
        glViewport(0, 0, width, height);

        // OpenGL Rendering related code
        glClearColor(0.2f, 0.3f, 0.3f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT);

				trans = glm::rotate(trans, (float)glfwGetTime(), glm::vec3(0.0f, 0.0f, 1.0f));
				transformLoc = glGetUniformLocation(ourShader.ID, "transform");
				glUniformMatrix4fv(transformLoc, 1, GL_FALSE, glm::value_ptr(trans));


        ourShader.use();

        glBindVertexArray(VAO);
        glDrawElements(GL_TRIANGLES, 6, GL_UNSIGNED_INT, 0);

        // Swap front and back buffers
        glfwSwapBuffers(window);

        // Poll for and process events
        glfwPollEvents();
    }

    glfwDestroyWindow(window);
    glfwTerminate();
    return 0;
}
