//
// Created by Göksu Güvendiren on 2019-05-14.
//

#include "Scene.hpp"


void Scene::buildBVH() {
    printf(" - Generating BVH...\n\n");
    this->bvh = new BVHAccel(objects, 1, BVHAccel::SplitMethod::NAIVE);
}

Intersection Scene::intersect(const Ray &ray) const
{
    return this->bvh->Intersect(ray);
}

void Scene::sampleLight(Intersection &pos, float &pdf) const
{
    float emit_area_sum = 0;
    for (uint32_t k = 0; k < objects.size(); ++k) {
        if (objects[k]->hasEmit()){
            emit_area_sum += objects[k]->getArea();
        }
    }
    float p = get_random_float() * emit_area_sum;
    emit_area_sum = 0;
    for (uint32_t k = 0; k < objects.size(); ++k) {
        if (objects[k]->hasEmit()){
            emit_area_sum += objects[k]->getArea();
            if (p <= emit_area_sum){
                objects[k]->Sample(pos, pdf);
                break;
            }
        }
    }
}

bool Scene::trace(
        const Ray &ray,
        const std::vector<Object*> &objects,
        float &tNear, uint32_t &index, Object **hitObject)
{
    *hitObject = nullptr;
    for (uint32_t k = 0; k < objects.size(); ++k) {
        float tNearK = kInfinity;
        uint32_t indexK;
        Vector2f uvK;
        if (objects[k]->intersect(ray, tNearK, indexK) && tNearK < tNear) {
            *hitObject = objects[k];
            tNear = tNearK;
            index = indexK;
        }
    }


    return (*hitObject != nullptr);
}

// Implementation of Path Tracing

Vector3f Scene::castRay(const Ray &ray, int depth) const
{
    // TO DO Implement Path Tracing Algorithm here
    
    Vector3f L_dir(0.0), L_indir(0.0);
    auto inter_ray = intersect(ray);
    if (!inter_ray.happened) {
        return this->backgroundColor;
    }

    if (inter_ray.obj->hasEmit()) {
        return inter_ray.m->getEmission();
    }

    auto wo = ray.direction;
    auto N = inter_ray.normal;
    auto p = inter_ray.coords;

    Intersection inter_light; 
    float pdf_light;
    sampleLight(inter_light, pdf_light);
    

    // This is wrong, because inter_light is not tagged "happened"
    const double eps = 5e-4;
    const double offset = 1e-3;
    auto x = inter_light.coords;
    auto ws = (x - p).normalized();

    auto NN = inter_light.normal;
    auto emit = inter_light.emit;
    auto inter_mid = intersect(Ray(p + N * offset, ws));
    
    bool hit_light = inter_mid.happened && inter_mid.distance - (x - p).norm() + offset + eps > 0.0;

    auto m = inter_ray.m;
    // assert(inter_mid.obj);
    if (hit_light) {
        auto d2 = dotProduct(x - p, x - p);
        L_dir = emit 
            * m->eval(wo, ws, N)
            * std::max(0.0f, dotProduct(ws, N))
            * std::max(0.0f, dotProduct(-ws, NN)) 
            / d2 
            / std::max(0.0001f, pdf_light);
    }

    if (get_random_float() <= RussianRoulette) {
        auto wi = (inter_ray.m->sample(wo, N)).normalized();
        auto inter_indir = intersect(Ray(p + N * offset, wi));
        if (inter_indir.happened && !inter_indir.obj->hasEmit()) {
            auto q = inter_indir.coords;
            auto Nq = inter_indir.normal;
            // bounding 
            L_indir = castRay(Ray(q + Nq * offset, wi), depth + 1) 
                * m->eval(wo, wi, N) 
                * std::max(0.0f, dotProduct(wi, N))
                / std::max(0.0001f, m->pdf(wo, wi, N)) 
                / RussianRoulette;
        }
    }
    return Vector3f::Max(Vector3f::Min(L_dir + L_indir, Vector3f(1.0f)), Vector3f(0.0f));
}
