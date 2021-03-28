//
//  Renderer.swift
//  Cornell Box
//
//  Created by Mark Lim Pak Mun on 09/02/2021.
//  Copyright © 2021 Mark Lim Pak Mun. All rights reserved.
//

import Cocoa
import MetalKit

class MetalView: MTKView, NSWindowDelegate {
    var computePipelineState: MTLComputePipelineState!
    var queue: MTLCommandQueue! = nil
    var u_time: Float = 0
    var u_timeBuffer: MTLBuffer!
    var u_mouseBuffer: MTLBuffer!

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
        let point = NSPoint(x: self.frame.width/2,
                            y: self.frame.height/2)
        setMouseCoords(point)
    }
    
    func buildPipeline() {
        // Load all the shader files with a metal file extension in the project
        guard let library = device!.newDefaultLibrary()
        else {
            fatalError("Could not load default library from main bundle")
        }

        // Create the render pipeline for the drawable render pass.
        let kernel = library.makeFunction(name: "compute")!
        
        do {
            computePipelineState = try self.device!.makeComputePipelineState(function: kernel)
        }
        catch {
            fatalError("Could not create compute pipeline state object: \(error)")
        }
    }

    func buildResources() {
        u_timeBuffer = device!.makeBuffer(length: MemoryLayout<Float>.stride,
                                          options: [])
        u_timeBuffer.label = "time uniform"
        u_mouseBuffer = device!.makeBuffer(length: MemoryLayout<float2>.stride,
                                           options: [])
        u_mouseBuffer.label = "mouse coords uniform"
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
        //print("draw")
        // Check the CAMetalDrawable instance for the current frame is renderable
        guard let drawable = currentDrawable
        else {
            return
        }
        updateTime()
        let commandBuffer = queue.makeCommandBuffer()
        let commandEncoder = commandBuffer.makeComputeCommandEncoder()

        commandEncoder.setComputePipelineState(computePipelineState)
        commandEncoder.setTexture(drawable.texture,
                                  at: 0)
        commandEncoder.setBuffer(u_timeBuffer,
                                 offset: 0,
                                 at: 0)
        commandEncoder.setBuffer(u_mouseBuffer,
                                 offset: 0,
                                 at: 1)

        let width = computePipelineState.threadExecutionWidth
        let height = computePipelineState.maxTotalThreadsPerThreadgroup / width
        let threadsPerThreadgroup = MTLSizeMake(width, height, 1)
        let threadgroupsPerGrid = MTLSizeMake((drawable.texture.width + width - 1) / width,
                                              (drawable.texture.height + height - 1) / height,
                                              1)
        // The method "dispatchThreadgroups:threadsPerThreadgroup" could also be used;
        // the call was introduced in macOS 10.13 (High Sierra)
        commandEncoder.dispatchThreadgroups(threadgroupsPerGrid,
                                            threadsPerThreadgroup: threadsPerThreadgroup)
        commandEncoder.endEncoding()
        commandBuffer.present(drawable)
        commandBuffer.commit()
    }
}
