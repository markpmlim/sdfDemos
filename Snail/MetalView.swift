//
//  Renderer.swift
//  Snail
//
//  Created by Mark Lim Pak Mun on 09/02/2021.
//  Copyright © 2021 Mark Lim Pak Mun. All rights reserved.
//

import Cocoa
import MetalKit

class MetalView: MTKView, NSWindowDelegate {
    var renderPipelineState: MTLRenderPipelineState!
    var queue: MTLCommandQueue! = nil
    var u_time: Float = 0
    var u_resolutionBuffer: MTLBuffer!
    var u_timeBuffer: MTLBuffer!
    var u_mouse =  float2(0,0)
    var u_mouseBuffer: MTLBuffer!

    var inputTexture1: MTLTexture!
    var inputTexture2: MTLTexture!
    var inputTexture3: MTLTexture!

    required init(coder: NSCoder) {
        super.init(coder: coder)
        guard let metalGPU = MTLCreateSystemDefaultDevice()
        else {
            fatalError("No Metal-aware GPU available")
        }
        self.device = metalGPU
        // Make the view's texture objects suitable for sampling and pixel r/w operations.
        self.framebufferOnly = false
        queue = self.device!.makeCommandQueue()
        buildPipeline()
        buildResources()
        inputTexture1 = loadTexture("Texture1")
        inputTexture2 = loadTexture("Texture2")
        inputTexture3 = loadTexture("Texture3")
    }
    
    func buildPipeline() {
        // Load all the shader files with a metal file extension in the project
        guard let library = device!.makeDefaultLibrary()
        else {
            fatalError("Could not load default library from main bundle")
        }

        // Create the render pipeline for the drawable render pass.
        let pipelineStateDescriptor = MTLRenderPipelineDescriptor()
        
        // Load the vertex program into the library
        let quadVertexProgram = library.makeFunction(name: "vertexShader")
        // Load the fragment program into the library
        let quadFragmentProgram = library.makeFunction(name: "fragmentShader")
        pipelineStateDescriptor.sampleCount = self.sampleCount
        pipelineStateDescriptor.vertexFunction = quadVertexProgram
        pipelineStateDescriptor.fragmentFunction = quadFragmentProgram
        pipelineStateDescriptor.colorAttachments[0].pixelFormat = self.colorPixelFormat

        do {
            renderPipelineState = try self.device!.makeRenderPipelineState(descriptor: pipelineStateDescriptor)
        }
        catch {
            fatalError("Could not create  render pipeline state object: \(error)")
        }
    }
    func buildResources() {
        u_resolutionBuffer = device!.makeBuffer(length: MemoryLayout<float2>.stride,
                                               options: [])
        u_timeBuffer = device!.makeBuffer(length: MemoryLayout<Float>.stride,
                                          options: [])
        u_timeBuffer.label = "time uniform"
        u_mouseBuffer = device!.makeBuffer(length: MemoryLayout<float2>.stride,
                                           options: [])
        u_mouseBuffer.label = "mouse coords uniform"
    }

    func loadTexture(_ fileName: String) -> MTLTexture {
        let options: [MTKTextureLoader.Option : Any] = [
            .textureUsage : NSNumber(value: MTLTextureUsage.shaderRead.rawValue),
            .textureStorageMode: NSNumber(value: MTLStorageMode.private.rawValue)
        ]
        let textureLoader = MTKTextureLoader(device: self.device!)
        var inputTexture: MTLTexture

        do {
            try inputTexture = textureLoader.newTexture(name: fileName,
                                                        scaleFactor: 1.0,
                                                        bundle: nil,
                                                        options: options)
        }
        catch {
            fatalError("Can't load the cube texture ")
        }
        return inputTexture
    }

    func updateTime() {
        u_time += Float(self.preferredFramesPerSecond)
        let bufferPointer = u_timeBuffer.contents()
        memcpy(bufferPointer, &u_time, MemoryLayout<Float>.stride)
    }

    func setMouseCoords(_ point: NSPoint) {
        var mouseLocation = float2(Float(point.x), Float(point.y))
        let bufferPointer = u_mouseBuffer.contents()
        memcpy(bufferPointer, &mouseLocation, MemoryLayout<float2>.stride)
    }

    override func mouseDown(with event: NSEvent) {
        let mousePoint = self.convert(event.locationInWindow, from: nil)
        setMouseCoords(mousePoint)
    }

    override func mouseDragged(with event: NSEvent) {
        let mousePoint = self.convert(event.locationInWindow, from: nil)
        setMouseCoords(mousePoint)
    }

    override func draw(_ dirtyRect: NSRect) {
        guard let drawable = currentDrawable
        else {
            return
        }
        guard let renderPassDescriptor = self.currentRenderPassDescriptor
        else {
            return
        }

        updateTime()
        guard let commandBuffer = queue.makeCommandBuffer()  else {
            return
        }
        guard let commandEncoder = commandBuffer.makeRenderCommandEncoder(descriptor:renderPassDescriptor) else {
            return
        }
        var resolution: float2 = float2(Float(self.drawableSize.width),
                                        Float(self.drawableSize.height))
        let bufferPointer = u_resolutionBuffer.contents()
        memcpy(bufferPointer, &resolution, MemoryLayout<float2>.stride)
        commandEncoder.setRenderPipelineState(renderPipelineState)
        commandEncoder.setFragmentTexture(inputTexture1,
                                  index: 0)
        commandEncoder.setFragmentTexture(inputTexture2,
                                  index: 1)
        commandEncoder.setFragmentTexture(inputTexture3,
                                  index: 2)
        commandEncoder.setFragmentBuffer(u_resolutionBuffer,
                                         offset: 0,
                                         index: 0)
        commandEncoder.setFragmentBuffer(u_timeBuffer,
                                         offset: 0,
                                         index: 1)
        commandEncoder.setFragmentBuffer(u_mouseBuffer,
                                         offset: 0,
                                         index: 2)
        commandEncoder.drawPrimitives(type: .triangleStrip,
                                      vertexStart: 0,
                                      vertexCount: 4)
        commandEncoder.endEncoding()
        commandBuffer.present(drawable)
        commandBuffer.commit()

    }
}
