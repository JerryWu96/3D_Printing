//Computational Fabrication Assignment #1
#include <iostream>
#include <vector>
#include "CompFab.h"
#include "Mesh.h"
#include <cfloat> // added DBL_MAX/DBL_MIN

//Triangle list (global)
typedef std::vector<CompFab::Triangle> TriangleList;

TriangleList g_triangleList;
CompFab::VoxelGrid *g_voxelGrid;

//Ray-Triangle Intersection
//Returns 1 if triangle and ray intersect, 0 otherwise
int rayTriangleIntersection(CompFab::Ray &ray, CompFab::Triangle &triangle)
{
    /********* ASSIGNMENT *********/
    /* Ray-Triangle intersection test: Return 1 if ray intersects triangle, 
     * 0 otherwise */

    // Step 1: compute normal vector of the current triangle.
    CompFab::Vec3 v1 = triangle.m_v1;
    CompFab::Vec3 v2 = triangle.m_v2;
    CompFab::Vec3 v3 = triangle.m_v3;
    CompFab::Vec3 n = (v1 - v2) % (v1 - v3); // normal vector
    
    // Step 2: Check if the ray and the triangle are parallel by using dot product.
    double denominator = ray.m_direction * n;
    if (denominator == 0.0) // perpendicular to the normal vector. infinite / 0 solutions
        return 0;
    // Step 3: check if t is valid
    double numerator = (triangle.m_v2 - ray.m_origin) * n;
    double t = numerator / denominator;
    if (t >= 0.0 && t < DBL_MAX) {
        // Step 4: check if the intersection point X is on the inside of all three edges
        CompFab::Vec3 X =
                ray.m_origin +
                CompFab::Vec3(t * ray.m_direction.m_x, t * ray.m_direction.m_y, t * ray.m_direction.m_z);
        CompFab::Vec3 e1 = v2 - v1; // three edges of the triangle
        CompFab::Vec3 e2 = v3 - v2;
        CompFab::Vec3 e3 = v1 - v3;
        CompFab::Vec3 x1 = X - v1;
        CompFab::Vec3 x2 = X - v2;
        CompFab::Vec3 x3 = X - v3;
        if ((e1 % x1 * n > 0) && (e2 % x2 * n > 0) && (e3 % x3 * n > 0))
            return 1;
    }
    return 0;

}

//Number of intersections with surface made by a ray originating at voxel and cast in direction.

// Added a parameter triangle.
int numSurfaceIntersections(CompFab::Vec3 &voxelPos, CompFab::Vec3 &dir)
{
    unsigned int numHits = 0;
    
/* Check and return the number of times a ray cast in direction dir,
* from voxel center voxelPos intersects the surface */

    CompFab::Ray ray(voxelPos, dir);//Vec3 &origin, Vec3 &direction
    for (auto& triangle : g_triangleList)
        if (rayTriangleIntersection(ray, triangle)) // 1: intersect!
            numHits++;

    return numHits;
}

bool loadMesh(char *filename, unsigned int dim)
{
    g_triangleList.clear();
    
    Mesh *tempMesh = new Mesh(filename, true);
    
    CompFab::Vec3 v1, v2, v3;

    //copy triangles to global list
    for(unsigned int tri =0; tri<tempMesh->t.size(); ++tri)
    {
        v1 = tempMesh->v[tempMesh->t[tri][0]];
        v2 = tempMesh->v[tempMesh->t[tri][1]];
        v3 = tempMesh->v[tempMesh->t[tri][2]];
        g_triangleList.push_back(CompFab::Triangle(v1,v2,v3));
    }

    //Create Voxel Grid
    CompFab::Vec3 bbMax, bbMin; // store min/max points of the entire Voxel Grid
    BBox(*tempMesh, bbMin, bbMax);
    
    //Build Voxel Grid
    double bbX = bbMax[0] - bbMin[0];
    double bbY = bbMax[1] - bbMin[1];
    double bbZ = bbMax[2] - bbMin[2];
    double spacing;
    
    if(bbX > bbY && bbX > bbZ)
    {
        spacing = bbX/(double)(dim-2);
    } else if(bbY > bbX && bbY > bbZ) {
        spacing = bbY/(double)(dim-2);
    } else {
        spacing = bbZ/(double)(dim-2);
    }
    
    CompFab::Vec3 hspacing(0.5*spacing, 0.5*spacing, 0.5*spacing);

    // define the voxelGrid with the min point and max spacing
    g_voxelGrid = new CompFab::VoxelGrid(bbMin-hspacing, dim, dim, dim, spacing);

    delete tempMesh;
    
    return true;
   
}

void saveVoxelsToObj(const char * outfile)
{
 
    Mesh box;
    Mesh mout;
    int nx = g_voxelGrid->m_dimX;
    int ny = g_voxelGrid->m_dimY;
    int nz = g_voxelGrid->m_dimZ;
    double spacing = g_voxelGrid->m_spacing;
    
    CompFab::Vec3 hspacing(0.5*spacing, 0.5*spacing, 0.5*spacing);
    
    for (int ii = 0; ii < nx; ii++) {
        for (int jj = 0; jj < ny; jj++) {
            for (int kk = 0; kk < nz; kk++) {
                if(!g_voxelGrid->isInside(ii,jj,kk)){
                        continue;
                }
                CompFab::Vec3 coord(0.5f + ((double)ii)*spacing, 0.5f + ((double)jj)*spacing, 0.5f+((double)kk)*spacing);
                CompFab::Vec3 box0 = coord - hspacing;
                CompFab::Vec3 box1 = coord + hspacing;
                makeCube(box, box0, box1);
                mout.append(box);
            }
        }
    }

    mout.save_obj(outfile);
}

int main()
{

    //dimension of voxel grid (e.g. 32x32x32)
    unsigned int dim = 32;
    std::cout << "Please type the directory of the object.\n" << "For example: ./sphere.obj\n";
    char* source;
    std::cin >> source;
    loadMesh(source, dim);

    CompFab::Vec3 voxelPos;
    CompFab::Vec3 direction(1.0,0.0,0.0);
    /*
     * Iterate over all voxels in g_voxelGrid and test whether they are inside our outside of the
     * surface defined by the triangles in g_triangleList
     * Cast ray, check if voxel is inside or outside
     * even number-> outside, while odd number -> inside
     */

    double spacing = g_voxelGrid -> m_spacing;
    for (int i = 0; i < g_voxelGrid -> m_dimX; i++)
        for (int j = 0; j < g_voxelGrid -> m_dimY; j++)
            for (int k = 0; k < g_voxelGrid -> m_dimZ; k++)
            {
                voxelPos.m_x = g_voxelGrid -> m_lowerLeft.m_x + ((double)i) * spacing;// + 0.5f;
                voxelPos.m_y = g_voxelGrid -> m_lowerLeft.m_y + ((double)j) * spacing;// + 0.5f;
                voxelPos.m_z = g_voxelGrid -> m_lowerLeft.m_z + ((double)k) * spacing;// + 0.5f;
                int numOfIntersect = numSurfaceIntersections(voxelPos, direction);
                if (numOfIntersect % 2 != 0) // odd, inside
                    g_voxelGrid -> m_insideArray[k*(g_voxelGrid -> m_dimX * g_voxelGrid -> m_dimY)+j* g_voxelGrid -> m_dimY + i] = true;
            }

    //Write  voxel data as obj
    char* output = "./output.obj";
    saveVoxelsToObj(output);
    delete g_voxelGrid;
}
