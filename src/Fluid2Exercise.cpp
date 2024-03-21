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

struct Emitter 
{
    Vector2 pos;
    Vector2 size;
    Vector2 velDir;
    float velMag;
    Vector3 color;
};

#define SMOKE_DISPENSER

#define DEMO_1
//#define DEMO_2
//#define DEMO_3


#define SUB_STEPS
//#define RUNGE_KUTTA_4

const uint NUM_SUBSTEPS = 10;
const double EPSILON = 1e-6;

#ifndef SMOKE_DISPENSER
#ifdef DEMO_1
const Emitter emitters[4] = {
    {Vector2(0.1, 0.1), Vector2(0.05f, 0.05f), Vector2(1.0f, 1.0f), 3.0f, Vector3(1.0f, 0.0f, 0.0f)},
    {Vector2(0.9, 0.9), Vector2(0.05f, 0.05f), Vector2(-1.0f, -1.0f), 3.0f, Vector3(0.0f, 1.0f, 0.0f)},
    {Vector2(0.9, 0.1), Vector2(0.05f, 0.05f), Vector2(-1.0f, 1.0f), 3.0f, Vector3(0.0f, 0.0f, 1.0f)},
    {Vector2(0.1, 0.9), Vector2(0.05f, 0.05f), Vector2(1.0f, -1.0f), 3.0f, Vector3(1.0f, 1.0f, 0.0f)}
};
#else // DEMO_1
#ifdef DEMO_2
const Emitter emitters[1] = {
    {Vector2(0.5, 0.1), Vector2(0.05f, 0.05f), Vector2(0.0f, 1.0f), 2.0f, Vector3(0.8f, 0.8f, 0.8f)}
};

#else // DEMO_2
#ifdef DEMO_3
const Emitter emitters[3] = {
    {Vector2(0.25, 0.1), Vector2(0.05f, 0.05f), Vector2(0.0f, 1.0f), 3.0f, Vector3(1.0f, 1.0f, 0.0f)},
    {Vector2(0.5, 0.1), Vector2(0.05f, 0.05f), Vector2(0.0f, 1.0f), 3.0f, Vector3(0.0f, 1.0f, 1.0f)},
    {Vector2(0.75, 0.1), Vector2(0.05f, 0.05f), Vector2(0.0f, 1.0f), 3.0f, Vector3(1.0f, 0.0f, 1.0f)}
};

#else // DEMO_3
const Emitter emitters[4] = {
    {Vector2(0.25, 0.1), Vector2(0.05f, 0.05f), Vector2(0.0, 1.0f), 2.0f, Vector3(1.0f, 0.0f, 0.0f)},
    {Vector2(0.75, 0.9), Vector2(0.05f, 0.05f), Vector2(0.0, -1.0f), 2.0f, Vector3(0.0f, 1.0f, 0.0f)},
    {Vector2(0.75, 0.1), Vector2(0.05f, 0.05f), Vector2(0.0, 1.0f), 2.0f, Vector3(0.0f, 0.0f, 1.0f)},
    {Vector2(0.25, 0.9), Vector2(0.05f, 0.05f), Vector2(0.0f, -1.0f), 2.0f, Vector3(1.0f, 1.0f, 0.0f)}};
#endif // DEMO_3
#endif // DEMO_2
#endif  // DEMO_1
#endif  // !SMOKE_DISPENSER

const Vector2 gDir = Vector2(0, 1);

double time = 0.0;

Vector3 colors[] = {
    Vector3(1.0f, 1.0f, 0.5f),
    Vector3(0.8f, 0.2f, 0.7f),
    Vector3(1.0f, 0.7f, 1.0f),
    Vector3(0.2f, 1.0f, 0.9f),
    Vector3(0.8f, 1.0f, 0.2f),
    Vector3(0.4f, 0.9f, 1.0f),
    Vector3(1.0f, 1.0f, 1.0f),
    Vector3(1.0f, 1.0f, 0.5f)  //(para volver al inicio)
};

Vector3 calculateTimeBasedColor(double time)
{
    int numColors = sizeof(colors) / sizeof(colors[0]);
    float cycleTime = 15.0f;  // Tiempo en segundos para completar un ciclo RGB completo
    float fraction = std::fmod(time / cycleTime, 1.0f) * (numColors - 1);
    int startIndex = static_cast<int>(fraction);
    int endIndex = (startIndex + 1) % numColors;
    float lerpFactor = fraction - startIndex;

    Vector3 StartColor = colors[startIndex];
    Vector3 EndColor = colors[endIndex];

    return Vector3(StartColor.x + (EndColor.x - StartColor.x) * lerpFactor,
                   StartColor.y + (EndColor.y - StartColor.y) * lerpFactor,
                   StartColor.z + (EndColor.z - StartColor.z) * lerpFactor);
}

template <typename T>
T biLerp(const Array2<T> &data, const Vector2 &posIJ)
{
    // Getting the 4 indexes
    const uint x0 = clamp(floor(posIJ.x), 0, data.getSize().x - 1);
    const uint x1 = clamp(x0 + 1, 0, data.getSize().x - 1);
    const uint y0 = clamp(floor(posIJ.y), 0, data.getSize().y - 1);
    const uint y1 = clamp(y0 + 1, 0, data.getSize().y - 1);

    // Setting the weights
    const float s = clamp(posIJ.x - static_cast<float>(x0), 0.0f, 1.0f);
    const float t = clamp(posIJ.y - static_cast<float>(y0), 0.0f, 1.0f);

    // getting the 4 values
    T const bottomLeft = data.getValue(x0, y0);
    T const bottomRight = data.getValue(x1, y0);
    T const topLeft = data.getValue(x0, y1);
    T const topRight = data.getValue(x1, y1);

    // returning the bilinear interpolation
    return bottomLeft * (1 - s) * (1 - t) + bottomRight * s * (1 - t) + topLeft * (1 - s) * t + topRight * s * t;
}

} // namespace
void Fluid2::fluidAdvection(const float dt)
{
    const float subtep = dt / NUM_SUBSTEPS;

    {
        Array2<Vector3> ink_copy = Array2<Vector3>(this->getInk());

        // Recorrer todas las celdas
        #pragma omp parallel for collapse(2)
        for (int i = 0; i < inkRGB.getSize().x; i++) {
            // For each column j
            for (int j = 0; j < inkRGB.getSize().y; j++) {
                // const auto x = grid get x(i, j)
                Vector2 pos = grid.getCellPos(Index2(i, j));

                float vx = (velocityX[Index2(i, j)] + velocityX[Index2(i + 1, j)]) * 0.5f;
                float vy = (velocityY[Index2(i, j)] + velocityY[Index2(i, j + 1)]) * 0.5f;

                Vector2 vel = Vector2(vx, vy);

#ifdef RUNGE_KUTTA_4

                Vector2 k1 = vel * dt;
                Vector2 pRK = pos - k1 * 0.5f;

                Vector2 k2 =
                    Vector2(biLerp(velocityX, grid.getCellIndex(pRK)), biLerp(velocityY, grid.getCellIndex(pRK))) * dt;
                pRK = pos - k2 * 0.5f;

                Vector2 k3 =
                    Vector2(biLerp(velocityX, grid.getCellIndex(pRK)), biLerp(velocityY, grid.getCellIndex(pRK))) * dt;
                pRK = pos - k3;

                Vector2 k4 =
                    Vector2(biLerp(velocityX, grid.getCellIndex(pRK)), biLerp(velocityY, grid.getCellIndex(pRK))) * dt;

                pos -= (k1 + 2 * k2 + 2 * k3 + k4) / 6.0f;

#else 
#ifdef SUB_STEPS
                for (int k = 0; k < NUM_SUBSTEPS; k++)
                {
	                pos -= subtep * vel;
                    
                        Vector2 ssIJ = grid.getCellIndex(pos);

                    float vx = biLerp(velocityX, ssIJ);
                    float vy = biLerp(velocityY, ssIJ);

                    vel = Vector2(vx, vy);
                }

#else // EULER
                pos -= dt * vel;
#endif // RUNGE_KUTTA_4 
#endif // SUB_STEPS
                const Vector2 ijprev = grid.getCellIndex(pos);

                Vector3 interpolated_ink = biLerp(ink_copy, ijprev);

                inkRGB[Index2(i, j)] = interpolated_ink;
            }
        }
    }

    {
        Array2<float> velX = Array2<float>(this->getVelocityX());
        Array2<float> velY = Array2<float>(this->getVelocityY());

        #pragma omp parallel for collapse(2)
        for (int i = 0; i < velocityX.getSize().x; i++) {
            for (int j = 0; j < velocityX.getSize().y; j++) {
                Vector2 pos = grid.getFacePosX(Index2(i, j));

                float vx = velX[Index2(i, j)];

                int i_y = clamp(i, 0, velocityY.getSize().x - 1);
                int j_y = j;
                int iminus1_y = clamp(i - 1, 0, velocityY.getSize().x - 1);
                int jminus1_y = clamp(j - 1, 0, velocityY.getSize().y - 1);

                float vy = bilerp(velY[Index2(i_y, j_y)], velY[Index2(iminus1_y, j_y)], velY[Index2(i_y, jminus1_y)], velY[Index2(iminus1_y, jminus1_y)], 0.5, 0.5);

                Vector2 vel = Vector2(vx, vy);

#ifdef RUNGE_KUTTA_4
                Vector2 k1 = vel * dt;
                Vector2 pRK = pos - k1 * 0.5f;

                Vector2 k2 = Vector2(biLerp(velX, grid.getFaceIndexX(pRK)), biLerp(velY, grid.getFaceIndexX(pRK))) *
                    dt;
                pRK = pos - k2 * 0.5f;

                Vector2 k3 = Vector2(biLerp(velX, grid.getFaceIndexX(pRK)), biLerp(velY, grid.getFaceIndexX(pRK))) *
                    dt;
                pRK = pos - k3;

                Vector2 k4 =
                    Vector2(biLerp(velocityX, grid.getFaceIndexX(pRK)), biLerp(velY, grid.getFaceIndexX(pRK))) *
                    dt;

                pos -= (k1 + 2 * k2 + 2 * k3 + k4) / 6.0f;

#else
#ifdef SUB_STEPS
                for (int k = 0; k < NUM_SUBSTEPS; k++)
                {
                    pos -= subtep * vel;

                    Vector2 ssIJ_x = grid.getFaceIndexX(pos);
                    Vector2 ssIJ_y = grid.getFaceIndexY(pos);

                    float vx = biLerp(velX, ssIJ_x);
                    float vy = biLerp(velY, ssIJ_y);

                    vel = Vector2(vx, vy);
                }
#else   // EULER
                pos -= dt * vel;
#endif
#endif  
                Vector2 ijprev = grid.getFaceIndexX(pos);


                float interpolated_vel = biLerp(velX, ijprev);

                // Update velocity
                velocityX[Index2(i, j)] = interpolated_vel;
            }
        }

        #pragma omp parallel for collapse(2)
        for (int i = 0; i < velocityY.getSize().x; i++) {
            for (int j = 0; j < velocityY.getSize().y; j++) {
                Vector2 pos = grid.getFacePosY(Index2(i, j));
                
                float vy =  velY[Index2(i,j)];

                int i_x = i;
                int j_x = clamp(j, 0, velocityX.getSize().y - 1);
                int iminus1_x = clamp(i - 1, 0, velocityX.getSize().x - 1);
                int jminus1_x = clamp(j - 1, 0, velocityX.getSize().y - 1);

                float vx = bilerp(velX[Index2(i_x, j_x)], velX[Index2(iminus1_x, j_x)], velX[Index2(i_x, jminus1_x)], velX[Index2(iminus1_x, jminus1_x)], 0.5, 0.5);
                
                Vector2 vel = Vector2(vx, vy);

                // RK - 4
#ifdef RUNGE_KUTTA_4
                Vector2 k1 = vel * dt;
                Vector2 pRK = pos - k1 * 0.5f;

                Vector2 k2 = Vector2(biLerp(velX, grid.getFaceIndexY(pRK)), biLerp(velY, grid.getFaceIndexY(pRK))) *
                    dt;
                pRK = pos - k2 * 0.5f;

                Vector2 k3 = Vector2(biLerp(velX, grid.getFaceIndexY(pRK)), biLerp(velY, grid.getFaceIndexY(pRK))) *
                    dt;
                pRK = pos - k3;

                Vector2 k4 = Vector2(biLerp(velX, grid.getFaceIndexY(pRK)), biLerp(velY, grid.getFaceIndexY(pRK))) *
                    dt;

                pos -= (k1 + 2 * k2 + 2 * k3 + k4) / 6.0f;
#else
#ifdef SUB_STEPS
                for (int k = 0; k < NUM_SUBSTEPS; k++)
				{
					pos -= subtep * vel;

					Vector2 ssIJ_x = grid.getFaceIndexX(pos);
                    Vector2 ssIJ_y = grid.getFaceIndexY(pos);

                    float vx = biLerp(velX, ssIJ_x);
                    float vy = biLerp(velY, ssIJ_y);

                    vel = Vector2(vx, vy);

					vel = Vector2(vx, vy);
				}
#else  // EULER
                pos -= dt * vel;
#endif
#endif  
                Vector2 ijprev = grid.getFaceIndexY(pos);

                float interpolated_vel = biLerp(velY, ijprev);

                velocityY[Index2(i, j)] = interpolated_vel;
            }
        }
    }
}

void Fluid2::fluidEmission()
{
    if (Scene::testcase >= Scene::SMOKE)
    {
#ifndef SMOKE_DISPENSER
        // Por cada uno de los emisores y fuerzo una velocidad concreta y una tinta.
        for (const auto& emitter : emitters)
        {
            #pragma omp parallel for collapse(2)
            for (int j = 0; j < inkRGB.getSize().y; ++j)
            {
                for (int i = 0; i < inkRGB.getSize().x; ++i)
                {
                    Index2 index = Index2(i, j);
                    const float xRel = static_cast<float>(i) / inkRGB.getSize().x;
                    const float yRel = static_cast<float>(j) / inkRGB.getSize().y;

                    if (std::abs(xRel - emitter.pos.x) < emitter.size.x / 2 &&
                        std::abs(yRel - emitter.pos.y) < emitter.size.y / 2)
                    {
						inkRGB[index] = emitter.color;
						velocityX[index] = emitter.velDir.x * emitter.velMag;
						velocityY[index] = emitter.velDir.y * emitter.velMag;
					}
				}
			}
		}   
 #else
        float pos = 0.5f;
        float size = 0.05f;
        float angle = std::fmod(time * 0.5f, 2 * M_PI);
        float anglePulse = std::fmod(time * 16, 2 * M_PI);

#ifdef DEMO_1
        float velmag = (std::sin(anglePulse) + 1) * 5;
        float velX = std::cos(angle) * velmag;
        float velY = std::sin(angle) * velmag;
#else
        float velmag = std::sin(anglePulse) * 10.0f;
        float velX = 1.0f * velmag;
        float velY = 0.0f;
#endif  // DEMO_1
        Vector3 color = calculateTimeBasedColor(time);

        #pragma omp parallel for collapse(2)
        for (int j = 0; j < inkRGB.getSize().y; ++j) {
            for (int i = 0; i < inkRGB.getSize().x; ++i) {
                Index2 index = Index2(i, j);
                const float xRel = static_cast<float>(i) / inkRGB.getSize().x;
                const float yRel = static_cast<float>(j) / inkRGB.getSize().y;

                float distanceToCenter = std::sqrt(std::pow(xRel - pos, 2) + std::pow(yRel - pos, 2));

                // Verifica si el punto está dentro del círculo
                if (distanceToCenter <= size) {
                    inkRGB[index] = color;
                    velocityX[index] = velX;
                    velocityY[index] = velY;
                }
            }
        }

        time += Scene::step;
#endif  //  SMOKE_DISPENSER

    }
}

void Fluid2::fluidVolumeForces(const float dt)
{
    
    if (Scene::testcase >= Scene::SMOKE) 
    {
        const float gravity = Scene::kGravity;

        #pragma omp parallel for collapse(2)
        for (int j = 0; j < grid.getSize().y; ++j) {
            for (int i = 0; i < grid.getSize().x; ++i) {
                Index2 index = Index2(i, j);
                Vector2 vel = Vector2(velocityX[index], velocityY[index]);

                // Gravity
                vel += gDir * gravity * dt;

                velocityX[index] = vel.x;
                velocityY[index] = vel.y;
            }
        }

        // Faltaría aplicar la fuerza en las posiciones velX[velX.x -1, j] y velY[i, velY.y -1]

        #pragma omp parallel for
        for (int j = 0; j < velocityX.getSize().y; ++j) {
            Index2 index = Index2(velocityX.getSize().x - 1, j);
            float vel = velocityX[index];

            vel = vel + gDir.x * gravity * dt;

            velocityX[index] = vel;
        }

        #pragma omp parallel for
        for (int i = 0; i < velocityY.getSize().x; ++i) {
            Index2 index = Index2(i, velocityY.getSize().y - 1);
            float vel = velocityY[index];

            vel = vel + gDir.y * gravity * dt;

            velocityY[index] = vel;
        }
    }
    
}

void Fluid2::fluidViscosity(const float dt)
{
    if (Scene::testcase >= Scene::SMOKE)
    {
        // Copy the velocity fields on temporary arrays
        const Array2<float> Vx_temp = Array2<float>(velocityX);
        const Array2<float> Vy_temp = Array2<float>(velocityY);

        // Precalculate some constants
        const float rho = Scene::kDensity;
        const float mu = Scene::kViscosity;
        const float dx = grid.getDx().x;
        const float dy = grid.getDx().y;
        const float dx2 = dx * dx;
        const float dy2 = dy * dy;
        const float mu_dt_over_rho = mu *(dt / rho);

        // Apply viscosity to velocityX
        {
            #pragma omp parallel for
            for (int i = 0; i < velocityX.getSize().x; i++) {
                for (int j = 0; j < velocityX.getSize().y; j++) {
                    uint iPlus1_x = min(i + 1, static_cast<int>(Vx_temp.getSize().x) - 1);
                    uint iMinus1_x = max(i - 1, 0);

                    uint jPlus1_y = min(j + 1, static_cast<int>(Vx_temp.getSize().y) - 1);
                    uint jMinus1_y = max(j - 1, 0);

					float xDerivative = (Vx_temp[Index2(iPlus1_x, j)] - 2 * Vx_temp[Index2(i, j)] + Vx_temp[Index2(iMinus1_x, j)]) / dx2;
                    float yDerivative = (Vx_temp[Index2(i, jPlus1_y)] - 2 * Vx_temp[Index2(i, j)] + Vx_temp[Index2(i, jMinus1_y)]) / dy2;

					velocityX[Index2(i, j)] += mu_dt_over_rho * (xDerivative + yDerivative);
                }
            }
        }

        // Apply viscosity to velocityY
        {
            #pragma omp parallel for
            for (int i = 0; i < velocityY.getSize().x; i++) {
                for (int j = 0; j < velocityY.getSize().y; j++) {
					uint iPlus1_x = min(i + 1, static_cast<int>(Vy_temp.getSize().x) - 1);
					uint iMinus1_x = max(i - 1, 0);
					uint jPlus1_y = min(j + 1, static_cast<int>(Vy_temp.getSize().y) - 1);
					uint jMinus1_y = max(j - 1, 0);

                    float xDerivative = (Vy_temp[Index2(iPlus1_x, j)] - 2 * Vy_temp[Index2(i, j)] + Vy_temp[Index2(iMinus1_x, j)]) / dx2;
                    float yDerivative = (Vy_temp[Index2(i, jPlus1_y)] - 2 * Vy_temp[Index2(i, j)] + Vy_temp[Index2(i, jMinus1_y)]) / dy2;

                    velocityY[Index2(i, j)] += mu_dt_over_rho * (xDerivative + yDerivative);
                }
            }
        }
    }
}

void Fluid2::fluidPressureProjection(const float dt)
{
    
    if (Scene::testcase >= Scene::SMOKE) 
    {
        // Precalculate some constants
        const float rho = Scene::kDensity;
        const float mu = Scene::kViscosity;
        const float dx = grid.getDx().x;
        const float dy = grid.getDx().y;
        const float dx2 = dx * dx;
        const float dy2 = dy * dy;

        const int nX = pressure.getSize().x;
        const int nY = pressure.getSize().y;
        const int size = nX * nY;

        // Vector b
        std::vector<double> b;
        b.resize(size);

        // Matriz A
        SparseMatrix<double> A(size, 5);

        // Vector de incognitas [del mismo tamaño]
        std::vector<double> p;
        p.resize(size, 0.0);
        
        // Forzar que todas las velocidades de 
        #pragma omp parallel for
        for (int i = 0; i < velocityX.getSize().y; ++i) 
        {
            velocityX[Index2(0, i)] = 0.0f;
            velocityX[Index2(velocityX.getSize().x - 1, i)] = 0.0f;
        }

        #pragma omp parallel for
        for (int j = 0; j < velocityY.getSize().x; ++j) 
        {
            velocityY[Index2(j, 0)] = 0.0f;
            velocityY[Index2(j, velocityY.getSize().y - 1)] = 0.0f;
        }

        // Las presiones estan en los centros de las celdas
        #pragma omp parallel for collapse(2)
        for (int j = 0; j < nY; ++j) 
        {
            for (int i = 0; i < nX; ++i)
            {
                // indice Lineal de A
                int idx = i + j * nX;

                // Rellenamos con la suma de las divergencias de las velocidades
                double u = (velocityX[Index2(i + 1, j)] - velocityX[Index2(i, j)]) / dx;
                double v = (velocityY[Index2(i, j + 1)] - velocityY[Index2(i, j)]) / dy;
                
                b[idx] = (-rho / dt) * (u + v);
            }
        }

        // Rellenamos la matriz A
        #pragma omp parallel for collapse(2)
        for (int j = 0; j < nY; ++j) {
            for (int i = 0; i < nX; ++i) {
                int idx = i + j * nX;

                double value = 2.0 / dx2 + 2.0 / dy2;

                // Aplicar modificaciones de valor para bordes
                value -= (i == 0 || i == nX - 1) / dx2;
                value -= (j == 0 || j == nY - 1) / dy2;

                // Conectar con nodos vecinos si no estamos en un borde
                if (i > 0)
                    A.add_to_element(idx, idx - 1, -1.0 / dx2);
                if (i < nX - 1)
                    A.add_to_element(idx, idx + 1, -1.0 / dx2);
                if (j > 0)
                    A.add_to_element(idx, idx - nX, -1.0 / dy2);
                if (j < nY - 1)
                    A.add_to_element(idx, idx + nX, -1.0 / dy2);

                // Establecer el valor del elemento diagonal principal + un pequeño epsilon para aumentar la estabilidad
                A.set_element(idx, idx, value + EPSILON);
            }
        }
        
        // Rellenamos el vector de presiones con el valor en el instante anterior   
        #pragma omp parallel for collapse(2)     
        for (int j = 0; j < nY; ++j) 
        {
            for (int i = 0; i < nX; ++i) 
            {
                p[i + j * nX] = pressure[Index2(i, j)];
            }
        }
       
        PCGSolver<double> PCG ;

        PCG.set_solver_parameters(1e-4, 100);
        

        double residual = 1e-4;
        int iterations = 100;
        PCG.solve(A, b, p, residual, iterations);

        // Rellenamos pressure con el resultado de montar el escena
        #pragma omp parallel for collapse(2)
        for (int j = 0; j < nY; ++j) {
            for (int i = 0; i < nX; ++i) {
                // indice Lineal de A
                int idx = i + j * nX;

                pressure[Index2(i, j)] = p[idx];
            }
        }
        // Aplicamos el gradiente de presión
        #pragma omp parallel for collapse(2)
        for (int j = 0; j < grid.getSize().y; ++j) {
            for (int i = 0; i < grid.getSize().x; ++i) {
				velocityX[Index2(i + 1, j)] -= (dt / rho) * ((pressure[Index2(i + 1, j)] - pressure[Index2(i, j)]) / dx);
				velocityY[Index2(i, j + 1)] -= (dt / rho) * ((pressure[Index2(i, j + 1)] - pressure[Index2(i, j)]) / dy);
			}
		}
    }
}


}  // namespace asa
