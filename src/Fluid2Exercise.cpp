#include "Scene.h"

#include "Numeric/PCGSolver.h"
#include "omp.h"

namespace asa
{
namespace
{
////////////////////////////////////////////////
// Add any reusable classes or functions HERE //
////////////////////////////////////////////////

#define NUM_SUBSTEPS 10

const Vector2 gDir = Vector2(0, 1);

 template <typename T>
T biLerp(const Array2<T> &colors, const Vector2 &posIJ)
{
    // Getting the 4 indexes
    const uint x0 = clamp(      floor(posIJ.x), 0, colors.getSize().x - 1);
    const uint x1 = clamp(              x0 + 1, 0, colors.getSize().x - 1);
    const uint y0 = clamp(      floor(posIJ.y), 0, colors.getSize().y - 1);
    const uint y1 = clamp(              y0 + 1, 0, colors.getSize().y - 1);

    // Setting the weights
    const float s = posIJ.x - x0;
    const float t = posIJ.y - y0;

    // getting the 4 values
    T const bottomLeft = colors.getValue(x0, y0);
    T const bottomRight = colors.getValue(x1, y0);
    T const topLeft = colors.getValue(x0, y1);
    T const topRight = colors.getValue(x1, y1);

    // returning the bilinear interpolation
    return bottomLeft * (1 - s) * (1 - t) + bottomRight * s * (1 - t) + topLeft * (1 - s) * t + topRight * s * t;
}
}  // namespace

// advection
void Fluid2::fluidAdvection(const float dt)
{
    {
        // Ink advecion HERE
        const Array2<Vector3> inkCopy = Array2<Vector3>(this->getInk());

        #pragma omp parallel for collapse(2) 
        for (int j = 0; j < inkRGB.getSize().y; ++j) 
        {
            for (int i = 0; i < inkRGB.getSize().x; ++i)
            {
                Index2 index = Index2(i, j);
                Vector2 pos = grid.getCellPos(index);
                Vector2  vel = Vector2(velocityX[index], velocityY[index]);
                // Euler
                // pos = pos - (vel * dt);
                
                // RK 4
                /*
                Vector2 k1 = vel * dt;
                Vector2 pRK = pos - k1 * 0.5f;

                Vector2 k2 =
                    Vector2(biLerp(velocityX, grid.getCellIndex(pRK)), biLerp(velocityY, grid.getCellIndex(pRK))) * dt;
                pRK= pos - k2 * 0.5f;

                Vector2 k3 =
                    Vector2(biLerp(velocityX, grid.getCellIndex(pRK)), biLerp(velocityY, grid.getCellIndex(pRK))) * dt;
                pRK = pos - k3;

                Vector2 k4 =
					Vector2(biLerp(velocityX, grid.getCellIndex(pRK)), biLerp(velocityY, grid.getCellIndex(pRK))) * dt;

                pRK = pos - (k1 + 2 * k2 + 2 * k3 + k4) / 6.0f;
                */


                // Substeps
                const float substep = dt / NUM_SUBSTEPS;
                for (uint k = 0; k < NUM_SUBSTEPS; ++k)
                {
					pos = pos - (vel * substep);
                    Vector2 ssIJ = grid.getCellIndex(pos);

                    float vx = biLerp(velocityX, ssIJ);
                    float vy = biLerp(velocityY, ssIJ);

                    vel = Vector2(vx, vy);
				}
                
				auto const prevIJ = grid.getCellIndex(pos);

				auto inkNext = biLerp(inkCopy, prevIJ);

				this->inkRGB[index] = inkNext;
            }
		}
    }

    {
        // Velocity acvection term HERE
    }
}



void Fluid2::fluidEmission()
{
    if (Scene::testcase >= Scene::SMOKE) {
        // Emitters contribution HERE
    }
}

void Fluid2::fluidVolumeForces(const float dt)
{
    if (Scene::testcase >= Scene::SMOKE) 
    {
        const float gravity = Scene::kGravity;

        #pragma omp parallel for collapse(2)
        for (int j = 0; j < velocityY.getSize().y; ++j)
        {
            for (int i = 0; i < velocityX.getSize().x; ++i)
            {
				Index2 index = Index2(i, j);
				Vector2 vel = Vector2(velocityX[index], velocityY[index]);

				// Gravity
                vel = vel + gDir * gravity;

				velocityX[index] = vel.x;
				velocityY[index] = vel.y;
			}
		}
    }
}

void Fluid2::fluidViscosity(const float dt)
{
    if (Scene::testcase >= Scene::SMOKE)
    {
        // Copy the velocity to a temporary grid
        const Array2<float> Vx = Array2<float>(velocityX);
        const Array2<float> Vy = Array2<float>(velocityY);

        // Precalculate some constants
        const float rho = Scene::kDensity;
        const float mu = Scene::kViscosity;
        const float dx = grid.getDx().x;
        const float dy = grid.getDx().y;
        const float dx2 = dx * dx;
        const float dy2 = dy * dy;

        #pragma omp parallel for collapse(2)
        for (int j = 0; j < Vx.getSize().y; ++j)
        {
            for (int i = 0; i < Vy.getSize().x; ++i)
            {
				Index2 index = Index2(i, j);
                Vector2 uij = Vector2(Vx[index], Vy[index]);

                Vector2 uiPlus1j = Vector2( Vx.getValue(Index2(i + 1, j)), Vy.getValue(Index2(i + 1, j)));
                Vector2 uiMinus1j = Vector2(Vx.getValue(Index2(i - 1, j)), Vy.getValue(Index2(i - 1, j)));
                Vector2 uijPlus1 = Vector2( Vx.getValue(Index2(i, j + 1)), Vy.getValue(Index2(i, j + 1)));
                Vector2 uijMinus1 = Vector2(Vx.getValue(Index2(i, j - 1)), Vy.getValue(Index2(i, j - 1)));

                Vector2 newUij = uij + (dt / rho) * (mu * ((uiPlus1j - (2 * uij) + uiMinus1j) / dx2 + (uijPlus1 - (2 * uij) + uijMinus1) / dy2));

                // Save the new velocity in the temporary grid
                 velocityX[index] = newUij.x;
                 velocityY[index] = newUij.y;
            }
        }
    }
}

void Fluid2::fluidPressureProjection(const float dt)
{
    if (Scene::testcase >= Scene::SMOKE) {
        // Incompressibility / Pressure term HERE
    }
}


}  // namespace asa
