# Name: Omar Madhani

import sys

import numpy
import math


# Constants assigned to the appropriate list index
# not to be modified

# Sphere
NAME = 0
RED = 1
GREEN = 2
BLUE = 3
KA = 4
KD = 5
KS = 6
KR = 7
N = 8

# Light
# NAME (same as sphere)
IR = 1
IG = 2
IB = 3

# (reference for function) inspired from provided ppm.cpp implementation 
def outputFile(outputFileName, resX, resY, backgroundRed, backgroundGreen, backgroundBlue, output):
    print(outputFileName)
    pixels: [float] = [0]*resX*resY*3
    file = open(outputFileName, "w")
    file.write("P3\n")
    file.write(f"{resX} {resY}\n")
    file.write("255\n")
    k = 0
    for i in range(resY):
        for j in range(resX):
            if output[k] == -1 and output[k+1] == -1 and output[k+2] == -1:
                pixels[k] = backgroundRed * 255
                pixels[k+1] = backgroundGreen * 255
                pixels[k+2] = backgroundBlue * 255
            else:
                pixels[k] = output[k] * 255
                pixels[k+1] = output[k+1] * 255
                pixels[k+2] = output[k+2] * 255  
            file.write(f"{pixels[k]} {pixels[k+1]} {pixels[k+2]}\n")
            k += 3
    file.close()

# used for ads lighting equation
def normalizeNegativeIntersectionPostion(intersectionPoint):
    negativeIntersectionPostion = -1*intersectionPoint
    negativeIntersectionPostion = numpy.delete(negativeIntersectionPostion, 3)
    normalizedNegativeIntersectionPostion = negativeIntersectionPostion / numpy.linalg.norm(negativeIntersectionPostion)
    return normalizedNegativeIntersectionPostion

# calculates the direction of the light
def lightDirection(intersectionPoint, light):
    lightRay = numpy.subtract(light, numpy.delete(intersectionPoint, 3))
    normalizedLightRay = lightRay / numpy.linalg.norm(lightRay)
    return normalizedLightRay
    
# calculates the surface normal at the intersection of a sphere
def surfaceNormal(intersectionPoint, sphereCenter, sphere):

    # we need to consider that the "sphere" is actually an ellipsoid
    difference = numpy.array([(intersectionPoint[0] - sphereCenter[0]) / (sphere[0][0])**2, (intersectionPoint[1] - sphereCenter[1]) / (sphere[1][1])**2, (intersectionPoint[2] - sphereCenter[2]) / (sphere[2][2])**2])
    difference = difference / numpy.linalg.norm(difference)
    return difference

# check if the light intersects with other spheres before the current sphere
def checkShadow(intersectionPoint, lightDirectionVector, currentSphere, sphereList, flippedNormal):
    for sphere in sphereList:
        if ((not (numpy.array_equal(currentSphere, sphere))) or flippedNormal):
            
            # inspired by CSC 305 Ray Tracing slide deck
            inverseModelMatrix = numpy.linalg.inv(sphere)
            origin = numpy.matmul(inverseModelMatrix, intersectionPoint)
            origin = numpy.delete(origin, 3)
            cValue = numpy.dot(origin, origin)-1
            newRay = numpy.matmul(inverseModelMatrix, numpy.array([lightDirectionVector[0], lightDirectionVector[1], lightDirectionVector[2], 0]))
            newRay = numpy.delete(newRay, 3)
            aValue = numpy.dot(newRay, newRay)
            bValue = 2*numpy.dot(newRay, origin)
            discriminantValue = (bValue**2) - (4 * aValue * cValue)
            if (discriminantValue >= 0):
                tValueAddition = (-1 * bValue + (math.sqrt(discriminantValue)))/(2*aValue)
                tValueSubtraction = (-1 * bValue - (math.sqrt(discriminantValue)))/(2*aValue)

                if min(tValueAddition, tValueSubtraction) > 0.000001:
                    return True
    return False

def calculateSpecularAndDiffuse(intersectionPoint, light, sphere, sphereList, flippedNormal, normal, lightIndex, lightAttributes, currentSphereRed, currentSphereGreen, currentSphereBlue, currentSphereKD, currentSphereKS, currentSphereSpecularExponent, output, outputIndex):

# only consider diffuse and specular lighting if the sphere is directly receiving light
    if not checkShadow(intersectionPoint, light, sphere, sphereList, flippedNormal):
        specularComponentRed = 0
        specularComponentGreen = 0
        specularComponentBlue = 0

        diffuseComponentRed = 0
        diffuseComponentGreen = 0
        diffuseComponentBlue = 0

        NdotL = max(0, numpy.dot(normal, light))

        # no need to consider diffuse if the coefficient is 0
        if not currentSphereKD == 0:
            diffuseComponentRed = currentSphereKD*lightAttributes[lightIndex][IR]*NdotL*currentSphereRed
            diffuseComponentGreen = currentSphereKD*lightAttributes[lightIndex][IG]*NdotL*currentSphereGreen
            diffuseComponentBlue = currentSphereKD*lightAttributes[lightIndex][IB]*NdotL*currentSphereBlue

        if NdotL > 0:
            
            # no need to consider specular if the coefficient is 0
            if not currentSphereKS == 0: 
                
                # inspired from previous WebGL assignment HTML code
                reflectValue = (((-1*light) - (2 * numpy.dot(normal, -1*light)*normal)))
                normalizedPosition = normalizeNegativeIntersectionPostion(intersectionPoint)
                RdotV = max(0, numpy.dot(reflectValue, normalizedPosition))
                specularComponentRed = currentSphereKS*lightAttributes[lightIndex][IR]*(RdotV**currentSphereSpecularExponent)
                specularComponentGreen = currentSphereKS*lightAttributes[lightIndex][IG]*(RdotV**currentSphereSpecularExponent)
                specularComponentBlue = currentSphereKS*lightAttributes[lightIndex][IB]*(RdotV**currentSphereSpecularExponent)
    
        output[outputIndex] += diffuseComponentRed + specularComponentRed
        output[outputIndex+1] += diffuseComponentGreen + specularComponentGreen
        output[outputIndex+2] += diffuseComponentBlue + specularComponentBlue

# need to send a ray to each pixel (resY x resX)
def ray(near, right, top, resX, resY, sphereList, sphereAttributes, lightMatrixList, lightAttributes, backgroundRed, backgroundGreen, backgroundBlue, ambientIR, ambientIG, ambientIB, outputFileName):
    output = [-1] * 600 * 600 * 3
    eye = numpy.array([0, 0, 0, 1])
    intersectionValues = [10000] * 600 * 600
    
    # attributes of spheres
    attributeIndex = 0
    
    for sphere in sphereList:

        inverseModelMatrix = numpy.linalg.inv(sphere)
        origin = numpy.matmul(inverseModelMatrix, eye)
        origin = numpy.delete(origin, 3)
        cValue = numpy.dot(origin, origin)-1

        currentSphereAttributes = sphereAttributes[attributeIndex]
        currentSphereKA = currentSphereAttributes[KA]
        currentSphereKD = currentSphereAttributes[KD]
        currentSphereKS = currentSphereAttributes[KS]
        currentSphereSpecularExponent = currentSphereAttributes[N]

        currentSphereRed = currentSphereAttributes[RED]
        currentSphereGreen = currentSphereAttributes[GREEN]
        currentSphereBlue = currentSphereAttributes[BLUE]

        currentSphereCenter = sphere[:, 3]
        
        tValueIndex = 0
        outputIndex = 0
        
        # we need to do calculations for each pixel
        for row in range(resY):
            for column in range(resX):
                
                # inspired by Ray Tracing CSC 305 slide deck

                ray = numpy.array([right*(((2*column)/resX) - 1), -1 * (top*(((2*row)/resY) - 1)), -1 * near, 0])
                
                newRay = numpy.matmul(inverseModelMatrix, ray)
                newRay = numpy.delete(newRay, 3)
                aValue = numpy.dot(newRay, newRay)
                bValue = 2*numpy.dot(newRay, origin)
        
                discriminantValue : float = (bValue **2) - (4 * aValue * cValue)
                
                # test if we have an intersection
                if discriminantValue >= 0:
                    tValueAddition = (-1 * bValue + (math.sqrt(discriminantValue)))/(2*aValue)
                    tValueSubtraction = (-1 * bValue - (math.sqrt(discriminantValue)))/(2*aValue)
                    
                    if tValueAddition > near and tValueSubtraction > near:
                        tValue = min(tValueSubtraction, tValueAddition)
                    elif tValueAddition > near or tValueSubtraction > near:
                        tValue = max(tValueSubtraction, tValueAddition)
                    else:
                        outputIndex += 3
                        continue
                    
                    # we only want to overwrite the current output pixel if the intersection point is closer to the eye
                    if tValue < intersectionValues[tValueIndex]:
                        intersectionValues[tValueIndex] = tValue

                        intersectionPoint = numpy.add(eye, (tValue)*ray)
                        
                        # intially set the output pixel to the ambient lighting of the current sphere, because external lights do not affect ambient lighting
                        output[outputIndex] = currentSphereKA*ambientIR*currentSphereRed
                        output[outputIndex+1] = currentSphereKA*ambientIG*currentSphereGreen
                        output[outputIndex+2] = currentSphereKA*ambientIB*currentSphereBlue
                        
                        normal = surfaceNormal(intersectionPoint, currentSphereCenter, sphere)

                        flippedNormal = False
                        if tValueAddition > near and tValueSubtraction < near:
                            normal = -1 * normal
                            flippedNormal = True
                       
                        
                        lightDirectionList = []
                        
                        # we need the direction of each light off of the current intersection point
                        for light in lightMatrixList:
                            lightDirectionList.append(lightDirection(intersectionPoint, light))
                            if flippedNormal == True:
                                break
                    
                        lightIndex = 0

                        # we need to sum up the diffuse and specular lighting contributions from each light
                        for l in lightDirectionList:
                            calculateSpecularAndDiffuse(intersectionPoint, l, sphere, sphereList, flippedNormal, normal, lightIndex, lightAttributes, currentSphereRed, currentSphereGreen, currentSphereBlue, currentSphereKD, currentSphereKS, currentSphereSpecularExponent, output, outputIndex)
                            lightIndex += 1

                outputIndex += 3
                tValueIndex += 1

        outputIndex = 0
        attributeIndex += 1

        # output at a pixel can be over 1, so we need to restrict each value of the array to 1
        clippedOutput = numpy.clip(output, -1, 1)
            
    outputFile(outputFileName, resX, resY, backgroundRed, backgroundGreen, backgroundBlue, clippedOutput)
    
# parses the input file
def readFile(fileName):
    near = None
    right = None
    top = None
    resX = None
    resY = None
    backgroundRed = None
    backgroundGreen = None
    backgroundBlue = None
    ambientIR = None
    ambientIG = None
    ambientIB = None
    outputFileName = None
    
    # contains the model matrix of each sphere
    sphereMatrixList = []
    sphereAttributes = [] 

    lightMatrixList = []
    lightAttributes = []


    file = open(fileName, "r")

    for line in file:
        content = line.strip().split()
        if len(content) == 0:
            continue
        elif content[0] == "NEAR":
            near = float(content[1])
        elif content[0] == "RIGHT":
            right = float(content[1])
        elif content[0] == "TOP":
            top = float(content[1])
        elif content[0] == "RES":
            resX = int(content[1])
            resY = int(content[2])
        elif content[0] == "SPHERE":
            sphereMatrix = numpy.array([[float(content[5]),0,0,float(content[2])],[0,float(content[6]),0,float(content[3])],[0,0,float(content[7]),float(content[4])],[0,0,0,1]])
            sphereMatrixList.append(sphereMatrix)
            sphereAttributes.append([content[1], float(content[8]), float(content[9]), float(content[10]), float(content[11]), float(content[12]), float(content[13]), float(content[14]), int(content[15])])
        elif content[0] == "LIGHT":
            lightMatrix = numpy.array([float(content[2]), float(content[3]), float(content[4])])
            lightMatrixList.append(lightMatrix)
            lightAttributes.append([content[1], float(content[5]), float(content[6]), float(content[7])])
        elif content[0] == "BACK":
            backgroundRed = float(content[1])
            backgroundGreen = float(content[2])
            backgroundBlue = float(content[3])
        elif content[0] == "AMBIENT":
            ambientIR = float(content[1])
            ambientIG = float(content[2])
            ambientIB = float(content[3])
        elif content[0] == "OUTPUT":
            outputFileName = content[1]

    file.close()
    ray(near, right, top, resX, resY, sphereMatrixList, sphereAttributes, lightMatrixList, lightAttributes, backgroundRed, backgroundGreen, backgroundBlue, ambientIR, ambientIG, ambientIB, outputFileName)


def main():
    readFile(sys.argv[1])

if __name__ == "__main__":
    main()
