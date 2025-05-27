import gmsh
import numpy as np

r = 1.2

gmsh.initialize()
gmsh.logger.start()
gmsh.model.add("test")
gmsh.option.setNumber("Geometry.Tolerance", 1e-9)
gmsh.option.setNumber("General.Terminal", 1)
gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
# Set debug
# gmsh.option.setNumber("General.Verbosity", 99)
gmsh.option.setNumber("Mesh.TransfiniteTri", 1)

gmsh.model.geo.synchronize()

p1 = gmsh.model.geo.addPoint(1, 0, 0, 10e-2)
p2 = gmsh.model.geo.addPoint(4, 0, 0, 10e-2)
p3 = gmsh.model.geo.addPoint(1, 1, 0, 10e-2)

l1 = gmsh.model.geo.addLine(p1, p2)
l2 = gmsh.model.geo.addLine(p2, p3)
l3 = gmsh.model.geo.addLine(p3, p1)

cl = gmsh.model.geo.addCurveLoop([l1, l2, l3])
s = gmsh.model.geo.addPlaneSurface([cl])

gmsh.model.geo.synchronize()


gmsh.option.setNumber("Mesh.RecombineAll", 1)
# gmsh.model.geo.synchronize()
# gmsh.model.mesh.setRecombine(2, s)
# gmsh.option.setNumber("Mesh.Algorithm", 8)  # 5 or 6
# gmsh.option.setNumber("Mesh.RecombinationAlgorithm", 2)  # 2 or 3

# print("Hmm0")
# gmsh.model.addPhysicalGroup(2, [s], 1)
# gmsh.model.setPhysicalName(2, 1, "surface")
# gmsh.model.geo.synchronize()
# gmsh.model.mesh.generate(2)
# gmsh.model.geo.synchronize()

n_angle = 3
angle = 10 * np.pi / 180  # 5 degrees
outDimTag = gmsh.model.geo.revolve(
    [(2, s)] , #, (1, l1), (1, l2), (1, l3)],
    0, 0, 0, 0, 1, 0, angle, numElements=[n_angle], recombine=True,
)
s1 = outDimTag[2][1]  # The first surface tag after revolution
s2 = outDimTag[3][1]  # The second surface tag after revolution
s3 = outDimTag[4][1]  # The third surface tag after revolution


gmsh.model.geo.synchronize()
volumes = gmsh.model.getEntities(3)
v = volumes[0][1]  # Get the first volume tag

gmsh.model.mesh.setTransfiniteCurve(l1, 15, "Progressiion", -r)
gmsh.model.mesh.setTransfiniteCurve(l2, 15, "Progressiion", r)
gmsh.model.mesh.setTransfiniteCurve(l3, 15, "Progressiion", -r)
gmsh.model.mesh.setTransfiniteSurface(s, "Right", [p3, p1, p2])

volumes = gmsh.model.getEntities(3)
gmsh.model.addPhysicalGroup(3, volumes[0], 1)
gmsh.model.setPhysicalName(3, 1, "volume")
gmsh.model.mesh.generate(3)
# gmsh.model.geo.synchronize()


# Need to correct the position of the internal nodes now.

allNodes = gmsh.model.mesh.getNodes(3, v, True, False)
internalNodes = gmsh.model.mesh.getNodes(3, v, False, False)

nodes1 = gmsh.model.mesh.getNodes(2, s1, True, False)
nodes2 = gmsh.model.mesh.getNodes(2, s2, True, False)
nodes3 = gmsh.model.mesh.getNodes(2, s3, True, False)

fixedNodes = list(nodes1[0]) + list(nodes2[0]) + list(nodes3[0])

elementTypes, elementTags, nodeTags = gmsh.model.mesh.getElements(3, v)
nElements = sum(len(et) for et in elementTags)


def filter_nodes(nodeTags, angle):
    """
    Returns the tags of nodes that only lie on the plane defined by the angle.
    """
    filteredNodeTags = []
    removedNodeTags = []
    for nTag in nodeTags:
        nodeCoord = gmsh.model.mesh.getNode(nTag)[0]
        nodeAngle = -np.arctan2(nodeCoord[2], nodeCoord[0])
        if np.isclose(nodeAngle, angle, atol=1e-6):
            filteredNodeTags.append(nTag)
        else:
            removedNodeTags.append(nTag)
    return filteredNodeTags, removedNodeTags


for iAngle in range(n_angle + 1):
    a = iAngle * angle / n_angle
    print(f"Angle: {np.rad2deg(a)} degrees")

    for _ in range(int(nElements/n_angle)):
        for i, eType in enumerate(elementTypes):
            for eTag in elementTags[i]:
                _, eNodes, _, _ = gmsh.model.mesh.getElement(eTag)
                eNodes, _ = filter_nodes(eNodes, a)
                if len(set(eNodes) & set(fixedNodes)) == 3 and eType == 5:
                    # print(f"  Updating element {eTag} with nodes {eNodes}")

                    eNodesSet = set(eNodes)

                    freeNode = list(eNodesSet - set(fixedNodes))[0]
                    cornerNodes = [n for n in eNodes if n != freeNode]

                    # Create the three vectors spanned by the three nodes
                    p1 = gmsh.model.mesh.getNode(cornerNodes[0])[0]
                    p2 = gmsh.model.mesh.getNode(cornerNodes[1])[0]
                    p3 = gmsh.model.mesh.getNode(cornerNodes[2])[0]

                    v1 = np.array(p2) - np.array(p1)
                    v2 = np.array(p3) - np.array(p1)
                    v3 = np.array(p3) - np.array(p2)

                    vecPairs = [(v1, v2, p1), (v1, v3, p2), (v2, v3, p3)]

                    for vec1, vec2, p in vecPairs:
                        if np.isclose(np.dot(vec1, vec2), 0, atol=1e-6):
                            newCoord = np.array(p) - vec1 + vec2
                            gmsh.model.mesh.setNode(freeNode, newCoord.tolist(), [])
                            fixedNodes.append(freeNode)
                            break
                        else:
                            continue
                        break
                    break



gmsh.write("transfinite_tri.msh")

gmsh.fltk.run()
gmsh.finalize()
