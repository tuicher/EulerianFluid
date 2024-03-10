#ifndef GRID_2_H
#define GRID_2_H

#include "Containers/Index2.h"
#include "Math/AABox2.h"

namespace asa
{
class Grid2
{
public:
    Grid2()
        : domain(Vector2(-1, -1), Vector2(1, 1))
        , size(100, 100)
    {
    }

    Grid2(const AABox2 &domain, const Index2 &size) { init(domain, size); }

    void init(const AABox2 &domain, const Index2 &size)
    {
        this->domain = domain;
        this->size = size;
    }

    Vector2 getDx() const
    {
        const Vector2 extent(domain.maxPosition - domain.minPosition);
        return Vector2(extent.x / float(size.x), extent.y / float(size.y));
    }

    AABox2 getDomain() const { return domain; }

    Index2 getSize() const { return size; }

    Index2 getSizeFaces(const uint axis) const {
        return axis ? getSizeFacesY() : getSizeFacesX();
    }

    Index2 getSizeFacesX() const { return Index2(size.x + 1, size.y); }

    Index2 getSizeFacesY() const { return Index2(size.x, size.y + 1); }

    Index2 getSizeNodes() const { return Index2(size.x + 1, size.y + 1); }

    Vector2 getCellPos(const Index2 &id) const
    {
        const Vector2 pos(domain.minPosition + Vector2(float(id.x) + 0.5f, float(id.y) + 0.5f) * getDx());
        return pos;
    }

    Vector2 getFacePos(const Index2 &id, const uint axis) const
    {
        return axis ? getFacePosY(id) : getFacePosX(id);
    }

    Vector2 getFacePosX(const Index2 &id) const
    {
        const Vector2 pos(domain.minPosition + Vector2(float(id.x), float(id.y) + 0.5f) * getDx());
        return pos;
    }

    Vector2 getFacePosY(const Index2 &id) const
    {
        const Vector2 pos(domain.minPosition + Vector2(float(id.x) + 0.5f, float(id.y)) * getDx());
        return pos;
    }

    Vector2 getNodePos(const Index2 &id) const
    {
        const Vector2 pos(domain.minPosition + Vector2(float(id.x), float(id.y)) * getDx());
        return pos;
    }

    Vector2 getCellIndex(const Vector2 &pos) const
    {
        const Vector2 ispos((pos - domain.minPosition) / getDx() - 0.5f);
        return ispos;
    }

    Vector2 getFaceIndex(const Vector2 &pos, const uint axis) const
    {
        return axis ? getFaceIndexY(pos) : getFaceIndexX(pos);
    }

    Vector2 getFaceIndexX(const Vector2 &pos) const
    {
        Vector2 ispos((pos - domain.minPosition) / getDx());
        ispos.y -= 0.5f;
        return ispos;
    }

    Vector2 getFaceIndexY(const Vector2 &pos) const
    {
        Vector2 ispos((pos - domain.minPosition) / getDx());
        ispos.x -= 0.5f;
        return ispos;
    }

private:
    AABox2 domain;
    Index2 size;
};
};  // namespace asa

#endif
