This project consists of several Raymarching demos posted on the website: www.shadertoy.com. The shaders had been converted to Apple's Metal Shading Language.


On the CPU side, code to setup a call to a Metal kernel function or vertex-fragment function pair is in the file MetalView.swift. The class MetalView is a sub-class of the class MTKView which in turn is a sub-class of NSView. So, the class inherits methods like mouseDown and mouseDragged from its ancestors.

The uniforms "u_mouse" and "u_time" are passed when a kernel function is called by the method draw. Another uniform "u_resolution" is passed when a vertex-fragment function pair. It's not necessary to pass this uniform because the size of the output texture of the kernel function can be determined.

When the display is rendered using a fragment function, a MTLRenderCommandEncoder object is created to pass the uniforms to the vertex and fragment functions. However, for this particular project none of the demos required parameters to be passed to the vertex function. The vertex function has imbedded data which is used to create a rectangular canvas in the range [-1.0, 1.0]. The output of the fragment function is written out to the currentDrawable's (an instance of CAMetalDrawable) texture.

For those demos where the display is rendered to an output texture by calling a (compute) kernel function, an instance of MTLComputeCommandEncoder is created. This object is used to passed the required uniforms to the kernel function. The programmer gets to decide how many pixels can be processed in parallel.
The output of the kernel function is written out to a texture.


General Notes on how to convert shaders from www.shadertoy.com to MSL


1) Whenever a vec3/float3 is passed as a parameter in a call to sample a texture, an environment map is required. Instead of loading six 2D graphics, we use a special kind of 2D graphic with dimensions 2:1 (aka equiRectangular map) to instantiate the environment map. A special function must be written to use access the 2D texture for environment lookups.

2) Functions like "normalize" can not be used on vectors in program scope because Metal considers such vectors to be in constant address space. Refer to the section "4) Address Spaces" of MSL specification manual. Here's a line from the manual:

The address space for a variable at program scope must be constant.

3) Matrices can not be initialized in program scope.

4) Because vectors, matrices etc at program scope are constants, variables must be passed as parameters either in the form of a user-specified struct or just a data type (e.g. float, float3) supported by MSL.

5) If a function returns a variable as an inout parameter, the "thread" keyword must be prepended to the name of the data type.




RaymarchingBasics:

Metal has a texture loader class named "MTKTextureLoader" which can load 2D and 3D graphics and instantiate an MTLTexture object.  In particular the following function:

	newTexture(name:scaleFactor:bundle:options:)

can be utilised to load equirectangular as well as six 2D cube maps. The graphics must be placed in a special folder named "Assets.xcassets". For this demo, we use the option "Texture Set" followed by a click-and-drag of 2D graphic into the empty space designated by XCode. The original shader uses a graphic with small dimensions. We searched and found a graphic of higher resolution which gives an identical background as the original. A good environment map is normally obtained by using graphics of resolution 2048:1024 and above.



sdPrimitives:

This is a Metal port of the source code for a reference implementation of basic primitives at the link:

https://www.shadertoy.com/view/Xds3zN

If a demo encounters the following error message when executing its kernel function,

Execution of the command buffer was aborted due to an error during execution. Invalid Resource (IOAF code 9)

the programmer should try to use a pair of vertex-fragment functions instead.

A kernel function ("compute shader" in OpenGL jargon) will process a number of pixels (specified by the programmer) in parallel whereas a fragment function will process one pixel at a time. As such pixel shaders have a highly restricted parallel programming model in which each fragment is computed completely independently of every other fragment. Resources can not be shared between each invocation of the fragment shader function.


A brief descripton of the 3D Signed Distance Field functions used to render these basic primitives is at the link:

https://www.iquilezles.org/www/articles/distfunctions/distfunctions.htm


The article is written by Inigo Quilez.




Snail: 

A shader written by Inigo Quilez. The original shader code was ported to Metal by writing a pair of vertex-fragment functions. 



CornellBox:

This demo tries to output a display to match a photograph of the rendered scene at the web-link

http://www.graphics.cornell.edu/online/box/box.jpg



Requirements:

XCode 9.x running under macOS 10.13.x.


Links:

www.shadertoy.com


