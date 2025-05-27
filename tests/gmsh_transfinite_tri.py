import gmsh
import numpy as np

r = 1.0

gmsh.initialize()
gmsh.logger.start()
gmsh.model.add("KSite49")
gmsh.option.setNumber("Geometry.Tolerance", 1e-9)
gmsh.option.setNumber("General.Terminal", 0)
gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
# Set debug
# gmsh.option.setNumber("General.Verbosity", 99)
gmsh.option.setNumber("Mesh.TransfiniteTri", 1)


p1 = gmsh.model.geo.addPoint(0, 0, 0, 10e-2)
p2 = gmsh.model.geo.addPoint(1, 0, 0, 10e-2)
p3 = gmsh.model.geo.addPoint(0, 3, 0, 10e-2)

l1 = gmsh.model.geo.addLine(p1, p2)
l2 = gmsh.model.geo.addLine(p2, p3)
l3 = gmsh.model.geo.addLine(p3, p1)

cl = gmsh.model.geo.addCurveLoop([l1, l2, l3])
s = gmsh.model.geo.addPlaneSurface([cl])

gmsh.model.geo.synchronize()

gmsh.model.mesh.setTransfiniteCurve(l1, 10, "Progressiion", -r)
gmsh.model.mesh.setTransfiniteCurve(l2, 10, "Progressiion", r)
gmsh.model.mesh.setTransfiniteCurve(l3, 10, "Progressiion", -r)
gmsh.model.mesh.setTransfiniteSurface(s, "Right", [p3, p1, p2])
# gmsh.option.setNumber("Mesh.RecombineAll", 1)
gmsh.model.mesh.setRecombine(2, s)

# print("Hmm0")
gmsh.model.addPhysicalGroup(2, [s], 1)
gmsh.model.setPhysicalName(2, 1, "surface")
gmsh.model.mesh.generate(2)
gmsh.model.geo.synchronize()

# print("Hmm1")
_ = gmsh.model.geo.revolve(
    [(2, s)], 0, 0, 0, 0, 1, 0, 5*np.pi/180, numElements=[1], recombine=True
)
gmsh.model.geo.synchronize()
volumes = gmsh.model.getEntities(3)
gmsh.model.addPhysicalGroup(3, volumes[0], 1)
gmsh.model.setPhysicalName(3, 1, "volume")
gmsh.model.geo.synchronize()
# print("Hmm2")
# gmsh.model.mesh.generate(3)
# print("Hmm3")


# Need to correct the position of the internal nodes now.

# allNodes = gmsh.model.mesh.getNodes(2, s, True, False)
# internalNodes = gmsh.model.mesh.getNodes(2, s, False, False)

# l1Nodes = gmsh.model.mesh.getNodes(1, l1, True, False)
# l2Nodes = gmsh.model.mesh.getNodes(1, l2, True, False)
# l3Nodes = gmsh.model.mesh.getNodes(1, l3, True, False)
# fixedNodes = list(l1Nodes[0]) + list(l3Nodes[0]) + list(l2Nodes[0])

# elementTypes, elementTags, nodeTags = gmsh.model.mesh.getElements(2, s)

# while len(fixedNodes) < len(allNodes[0]):
# # for _ in range(10):
#     for i, eType in enumerate(elementTypes):
#         for eTag in elementTags[i]:
#             _, eNodes, _, _ = gmsh.model.mesh.getElement(eTag)

#             # Check if three of the elements nodes are in fixedNodes:
#             if len(set(eNodes) & set(fixedNodes)) == 3 and eType == 3:

#                 eNodesSet = set(eNodes)
#                 freeNode = list(eNodesSet - set(fixedNodes))[0]
#                 cornerNodes = [n for n in eNodes if n != freeNode]

#                 # Create the three vectors spanned by the three nodes
#                 p1 = gmsh.model.mesh.getNode(cornerNodes[0])[0]
#                 p2 = gmsh.model.mesh.getNode(cornerNodes[1])[0]
#                 p3 = gmsh.model.mesh.getNode(cornerNodes[2])[0]

#                 v1 = np.array(p2) - np.array(p1)
#                 v2 = np.array(p3) - np.array(p1)
#                 v3 = np.array(p3) - np.array(p2)

#                 # Find the two that are closest to being orthogonal,
#                 # i.e we need to check three pairs

#                 vecPairs = [(v1, v2), (v1, v3), (v2, v3)]

#                 for vec1, vec2 in vecPairs:
#                     if np.isclose(np.dot(vec1, vec2), 0, atol=1e-6):
#                         # We can use these two vectors to calculate the new position
#                         # of the free node
#                         newCoord = np.array(p1) + vec1 + vec2
#                         gmsh.model.mesh.setNode(freeNode, newCoord.tolist(), [])
#                         fixedNodes.append(freeNode)
#                         # gmsh.model.geo.synchronize()
#                         break



                # Now we can calculate the new position
                # pNew = (cloesestNodes[1] - closestNodes[-1])
                # p0 = gmsh.model.mesh.getNode(closestNodes[0][0])[0]
                # p1 = gmsh.model.mesh.getNode(closestNodes[1][0])[0]
                # p2 = gmsh.model.mesh.getNode(closestNodes[2][0])[0]

                # newCoord = p0 + (p1 - p2)
                # gmsh.model.mesh.setNode(freeNode, newCoord, [])
                # fixedNodes.append(freeNode)
                # gmsh.model.geo.synchronize()
                # if eTag not in range(40, 48):
                #     gmsh.model.mesh.setVisibility([eTag], 1)
            
                #     break


# for eType, eTag in zip(elementTypes, elementTags):
#     print(f"Element type: {eType}, Element tag: {eTag}")
#     eNodes = gmsh.model.mesh.getElementEdgeNodes(eType, eTag)
#     print("Element tag:", eTag, "Nodes:", eNodes)

# Revolve the mesh to make it 3D:


# Export to .msh file

# gmsh.write("transfinite_tri.msh")

gmsh.fltk.run()
gmsh.finalize()
