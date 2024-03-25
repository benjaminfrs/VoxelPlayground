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
#include <glm/gtx/string_cast.hpp>
#include "camera.h"
#include "block.hpp"
#include "chunk.hpp"
#include <simplex/SimplexNoise.cpp>
#include <noise/noise.h>
#include <noise/noiseutils.h>
//#include "noiseutils.cpp"

void processInput(GLFWwindow *window);
void framebuffer_size_callback(GLFWwindow* window, int width, int height);
void mouse_callback(GLFWwindow* window, double xpos, double ypos);
void scroll_callback(GLFWwindow* window, double xoffset, double yoffset);

// settings
static const unsigned int SCR_WIDTH = 800;
static const unsigned int SCR_HEIGHT = 600;
bool CURSOR_MOVING = false;

Camera cam(glm::vec3(0, 68, 68));
//Time globals...
float deltaTime = 0.0f;
float lastFrame = 0.0f;

float lastX = SCR_WIDTH / 2.0f;
float lastY = SCR_HEIGHT / 2.0f;
bool firstMouse = true;

//GLFWwindow* initializeOpenGL(){
//    glfwInit();
//    GLFWwindow* window;
//
//    // Without these two hints, nothing above OpenGL version 2.1 is supported
//    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GLFW_TRUE);
//    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
//    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
//    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
//
//    // Create a windowed mode window and its OpenGL context
//    window = glfwCreateWindow(SCR_WIDTH, SCR_HEIGHT, "Running OpenGL on Mac", NULL, NULL);
//    if (!window)
//    {
//        std::cout << "Window failed to create";
//        glfwTerminate();
//        return -1;
//    }
//    // Make the window's context current
//    glfwMakeContextCurrent(window);
//    glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);
//    glfwSetCursorPosCallback(window, mouse_callback);
//    glfwSetScrollCallback(window, scroll_callback);
//    //------------------------------------Init and Configure
//    //
//		// Dark blue background
//    glClearColor(0.0f, 0.0f, 0.4f, 0.0f);
//    glEnable(GL_DEPTH_TEST);
//    glDepthFunc(GL_LESS);
//    glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);
//
//
//    return window;
//
//}
int main(void)
{
    //std::cout << "81" << std::endl;
    //GLFWwindow* window;
    glfwInit();
    // Without these two hints, nothing above OpenGL version 2.1 is supported
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GLFW_TRUE);

    // Create a windowed mode window and its OpenGL context
    GLFWwindow* window = glfwCreateWindow(SCR_WIDTH, SCR_HEIGHT, "Running OpenGL on Mac", NULL, NULL);
    if (window == NULL)
    {
        std::cout << "Window failed to create";
        glfwTerminate();
        return -1;
    }
    //std::cout << window;
    // Make the window's context current
    //int *width = nullptr;
    //int *height = nullptr;
    //glfwGetWindowSize(window, width, height);
    //std::cout << *width << "x" << *height << std::endl;
    //std::cout << width << "||" << height << std::endl;
    glfwMakeContextCurrent(window);
    glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);
    glfwSetCursorPosCallback(window, mouse_callback);
    glfwSetScrollCallback(window, scroll_callback);
    //------------------------------------Init and Configure
    //
		// Dark blue background
    glClearColor(0.0f, 0.0f, 0.4f, 0.0f);
    glEnable(GL_DEPTH_TEST);
    glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);

    //window = initializeOpenGL();

    //Build and compile shaders
    Shader ourShader("TransformVertexShader.vertexshader", "ColorFragmentShader.fragmentshader"); 
    Shader lightingShader("LightingVertexShader.vertexshader", "LightingFragmentShader.fragmentshader"); 
    

    //lighting configuration
		// lighting
		//glm::vec3 lightPos(1.0f, 66.0f, 66.0f);

		//glm::vec3 lightColor(0.33f, 0.42f, 0.18f);
		glm::vec3 objectColor(0.1f, 0.5f, 0.1f);
		//glm::vec3 result = lightColor * toyColor; // = (0.33f, 0.21f, 0.06f);
		//float ambientStrength = 0.1;
		//glm::vec3 ambient = ambientStrength * lightColor;
		//end lighting config
		
    unsigned int VBO, VAO;
    glGenVertexArrays(1, &VAO);
    glGenBuffers(1, &VBO);

    glBindVertexArray(VAO);

    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glBufferData(GL_ARRAY_BUFFER, sizeof(CUBE::vertices), CUBE::vertices, GL_STATIC_DRAW);

    // position attribute
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);
    // normal attribute
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)(3 * sizeof(float)));
    glEnableVertexAttribArray(1);


		unsigned int lightVAO;
		glGenVertexArrays(1, &lightVAO);
		glBindVertexArray(lightVAO);
		// we only need to bind to the VBO, the container's VBO's data already contains the data.
		glBindBuffer(GL_ARRAY_BUFFER, VBO);
		// set the vertex attribute 
		glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)0);
		glEnableVertexAttribArray(0);

    //glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
		//glEnable(GL_CULL_FACE);
		//glCullFace(GL_BACK);


    ChunkList chunks;
		//std::cout << chunk;
    //std::cout << chunk;
    //std::array<std::array<int, chunk.CHUNK_SIZE>, chunk.CHUNK_SIZE> heightMap{{0}};
    ////std::cout << "testttttt" << heightMap[10][10];
    //for(float i=0.0f;i<chunk.CHUNK_SIZE-1;i++){
    //  for(float j=0.0f;j<chunk.CHUNK_SIZE-1;j++){
    //    heightMap[int(i)][int(j)] = int(SimplexNoise::noise(i, j)*chunk.CHUNK_SIZE/2 + chunk.CHUNK_SIZE/2);
    //    //std::cout << "height: " << heightMap[int(i)][int(j)] << std::endl;
    //}}
    module::Perlin myModule;
    utils::NoiseMap heightMapTest;
		utils::NoiseMapBuilderPlane heightMapTestBuilder;
    std::cout << "133" << std::endl;
    double lowerAngleBound{0.0}, upperAngleBound {4.0}, lowerHeightBound {0.0}, upperHeightBound {4.0};
		heightMapTestBuilder.SetBounds (lowerAngleBound, upperAngleBound, lowerHeightBound, upperHeightBound);
		heightMapTestBuilder.SetDestSize (256, 256);
		heightMapTestBuilder.SetDestNoiseMap (heightMapTest);
		heightMapTestBuilder.SetSourceModule (myModule);
    std::cout << "133" << std::endl;
		heightMapTestBuilder.Build();
		//for(int i=0;i<chunk.CHUNK_SIZE;i++){
		//	for(int j=0;j<chunk.CHUNK_SIZE;j++){
		//	  //std::cout << i << j << std::endl;
		//	  //std::cout << heightMapTest.GetValue(i, j);
		//}}

    chunks.setChunkVisibility(heightMapTest);
    //chunks.debugVis();
		//chunk.debugMakeChunkVisible();
		//chunk.findSurroundedBlocks();
    lowerAngleBound += 4.0;
    upperAngleBound += 4.0;
    lowerHeightBound += 4.0;
    upperHeightBound += 4.0;
    heightMapTestBuilder.SetBounds(lowerAngleBound, upperAngleBound, lowerHeightBound, upperHeightBound);
    heightMapTestBuilder.Build();

		glm::vec3 lightPos(0.0f, 33, 33);
    while (!glfwWindowShouldClose(window))
    {
        float currentFrame = static_cast<float>(glfwGetTime());
        deltaTime = currentFrame - lastFrame;
        lastFrame = currentFrame;

        processInput(window);
				lightPos.x = 1.0f + sin(glfwGetTime()) * 2.0f;
				lightPos.z = sin(glfwGetTime() / 2.0f) * 1.0f;

        glClearColor(0.1f, 0.1f, 0.1f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        ourShader.use();
        //ourShader.setVec3("lightPos", cam.Position-glm::vec3(0.0f, 0.0f, -2.0f));
        ourShader.setVec3("lightPos", lightPos); 
        ourShader.setVec3("viewPos", cam.Position);
        ourShader.setVec3("objectColor", objectColor);
        ourShader.setVec3("lightColor",  1.0f, 1.0f, 1.0f);

        ourShader.setMat4("projection", glm::perspective(glm::radians(cam.Zoom), (float)SCR_WIDTH / (float)SCR_HEIGHT, 0.1f, 100.0f));
        ourShader.setMat4("view", cam.GetViewMatrix());


        glBindVertexArray(VAO);
        chunks.renderChunks(ourShader, glm::mat4 {1.0f});
        //ourShader.setMat4("model", glm::mat4(1.0f));
        //glDrawArrays(GL_TRIANGLES, 0, 36);

        // also draw the lamp object
        lightingShader.use();
        lightingShader.setMat4("projection", glm::perspective(glm::radians(cam.Zoom), (float)SCR_WIDTH / (float)SCR_HEIGHT, 0.1f, 100.0f));
        lightingShader.setMat4("view", cam.GetViewMatrix());
        glm::mat4 model{1.0f};
        model = glm::translate(model, lightPos);
        //model = glm::scale(model, glm::vec3(0.2f)); // a smaller cube
        lightingShader.setMat4("model", model);

        glBindVertexArray(lightVAO);
        glDrawArrays(GL_TRIANGLES, 0, 36);

				//if(!CURSOR_MOVING) {glfwSetCursorPos(window, SCR_WIDTH/2, SCR_HEIGHT/2);}

        glfwSwapBuffers(window);
        glfwPollEvents();
    }

    //std::cout << 0x101010101000000  << std::endl;
    //std::cout << window << std::endl;
    //std::cout << std::endl << "terminating";
    glfwDestroyWindow(window);
    //std::cout << std::endl << "cleanedup";
    std::cout << "before term";
    glfwTerminate();
    //std::cout << std::endl << "cleanedup";
    return 0;
}


void processInput(GLFWwindow *window)
{
    if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS)
        glfwSetWindowShouldClose(window, true);

    if (glfwGetKey(window, GLFW_KEY_W) == GLFW_PRESS)
        cam.ProcessKeyboard(FORWARD, deltaTime);
    if (glfwGetKey(window, GLFW_KEY_S) == GLFW_PRESS)
        cam.ProcessKeyboard(BACKWARD, deltaTime);
    if (glfwGetKey(window, GLFW_KEY_A) == GLFW_PRESS)
        cam.ProcessKeyboard(LEFT, deltaTime);
    if (glfwGetKey(window, GLFW_KEY_D) == GLFW_PRESS)
        cam.ProcessKeyboard(RIGHT, deltaTime);
}

void mouse_callback(GLFWwindow* window, double xposIn, double yposIn)
{
    float xpos = static_cast<float>(xposIn);
    float ypos = static_cast<float>(yposIn);

    if (firstMouse)
    {
        lastX = xpos;
        lastY = ypos;
        firstMouse = false;
    }

    float xoffset = xpos - lastX;
    float yoffset = lastY - ypos; // reversed since y-coordinates go from bottom to top

    lastX = xpos;
    lastY = ypos;

    cam.ProcessMouseMovement(xoffset, yoffset);
}

void scroll_callback(GLFWwindow* window, double xoffset, double yoffset)
{
    cam.ProcessMouseScroll(static_cast<float>(yoffset));
}

void framebuffer_size_callback(GLFWwindow* window, int width, int height)
{
    // make sure the viewport matches the new window dimensions; note that width and 
    // height will be significantly larger than specified on retina displays.
    glViewport(0, 0, width, height);
}
