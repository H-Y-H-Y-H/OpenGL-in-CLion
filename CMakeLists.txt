cmake_minimum_required(VERSION 3.12)
project(OpenGL)
set(CMAKE_CXX_STANDARD 14)
set(GLEW_H /usr/local/Cellar/glew/2.1.0/include/GL)
set(GLFW_H /usr/local/Cellar/glfw/3.3/include/GLFW)

#I put the included header document in my project dictionary directly.
set(GLAD_H /Users/yuhang/CLionProjects/EA_Project/include/glad)
set(KH_H /Users/yuhang/CLionProjects/EA_Project/include/KHR)
include_directories(${GLEW_H} ${GLFW_H} ${GLAD_H} ${KH_H})

# Add target link
set(GLEW_LINK /usr/local/Cellar/glew/2.1.0/lib/libGLEW.2.1.dylib)
set(GLFW_LINK /usr/local/Cellar/glfw/3.3/lib/libglfw.3.dylib)
link_libraries(${OPENGL} ${GLEW_LINK} ${GLFW_LINK})

# Execute compile command
# pay attention to the source path
set(SOURCE_FILES "src/glad.c" "main.cpp")
add_executable(OpenGL ${SOURCE_FILES})

# mac下这步很重要
if (APPLE)
    target_link_libraries(OpenGL "-framework OpenGL")
    target_link_libraries(OpenGL "-framework GLUT")
endif()
