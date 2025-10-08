#include <algorithm>
#include <cassert>
#include "BVH.hpp"

BVHAccel::BVHAccel(std::vector<Object*> p, int maxPrimsInNode,
                   SplitMethod splitMethod)
    : maxPrimsInNode(std::min(255, maxPrimsInNode)), splitMethod(splitMethod),
      primitives(std::move(p))
{
    time_t start, stop;
    time(&start);
    if (primitives.empty())
        return;

    root = recursiveBuild(primitives);

    time(&stop);
    double diff = difftime(stop, start);
    int hrs = (int)diff / 3600;
    int mins = ((int)diff / 60) - (hrs * 60);
    int secs = (int)diff - (hrs * 3600) - (mins * 60);

    printf(
        "\rBVH Generation complete: \nTime Taken: %i hrs, %i mins, %i secs\n\n",
        hrs, mins, secs);
}



BVHBuildNode* BVHAccel::recursiveBuild(std::vector<Object*> objects)
{
    BVHBuildNode* node = new BVHBuildNode();

    // Compute bounds of all primitives in BVH node
    Bounds3 bounds;
    for (int i = 0; i < objects.size(); ++i)
        bounds = Union(bounds, objects[i]->getBounds());
    switch (splitMethod) {
        case SplitMethod::NAIVE:

        if (objects.size() == 1) {
            // Create leaf _BVHBuildNode_
            node->bounds = objects[0]->getBounds();
            node->objects = {objects[0]};
            node->left = nullptr;
            node->right = nullptr;
            return node;
        }
        else if (objects.size() == 2) {
            node->left = recursiveBuild(std::vector{objects[0]});
            node->right = recursiveBuild(std::vector{objects[1]});

            node->bounds = Union(node->left->bounds, node->right->bounds);
            return node;
        }
        else {
            Bounds3 centroidBounds;
            for (int i = 0; i < objects.size(); ++i)
                centroidBounds =
                    Union(centroidBounds, objects[i]->getBounds().Centroid());
            int dim = centroidBounds.maxExtent();
            switch (dim) {
            case 0:
                std::sort(objects.begin(), objects.end(), [](auto f1, auto f2) {
                    return f1->getBounds().Centroid().x <
                        f2->getBounds().Centroid().x;
                });
                break;
            case 1:
                std::sort(objects.begin(), objects.end(), [](auto f1, auto f2) {
                    return f1->getBounds().Centroid().y <
                        f2->getBounds().Centroid().y;
                });
                break;
            case 2:
                std::sort(objects.begin(), objects.end(), [](auto f1, auto f2) {
                    return f1->getBounds().Centroid().z <
                        f2->getBounds().Centroid().z;
                });
                break;
            }

            auto beginning = objects.begin();
            auto middling = objects.begin() + (objects.size() / 2);
            auto ending = objects.end();

            auto leftshapes = std::vector<Object*>(beginning, middling);
            auto rightshapes = std::vector<Object*>(middling, ending);

            assert(objects.size() == (leftshapes.size() + rightshapes.size()));

            node->left = recursiveBuild(leftshapes);
            node->right = recursiveBuild(rightshapes);

            node->bounds = Union(node->left->bounds, node->right->bounds);
        }

        break;
        case SplitMethod::SAH:

        if (objects.size() == 1) {
            // Create leaf _BVHBuildNode_
            node->bounds = objects[0]->getBounds();
            node->objects = {objects[0]};
            node->left = nullptr;
            node->right = nullptr;
            return node;
        }
        else if (objects.size() == 2) {
            node->left = recursiveBuild(std::vector{objects[0]});
            node->right = recursiveBuild(std::vector{objects[1]});

            node->bounds = Union(node->left->bounds, node->right->bounds);
            return node;
        }
        else {
            Bounds3 centroidBounds;
            for (int i = 0; i < objects.size(); ++i)
                centroidBounds =
                    Union(centroidBounds, objects[i]->getBounds().Centroid());
            
            const int SAH_BUCKETS = 16;
            const double trav_cost = 0.125;

            int total_n = objects.size();
            double inv_total_area = 1.0 / bounds.SurfaceArea();

            Bounds3 prefix[SAH_BUCKETS];
            Bounds3 suffix[SAH_BUCKETS];

            int buckets[SAH_BUCKETS];

            int partion_dim = -1;
            int partion_bucket = -1;
            double inv_axis_extend;
            double buckets_step;
            double sah_cost = std::numeric_limits<double>::max();

            for (int dim = 0; dim < 3; dim++) {
                inv_axis_extend = 1.0 / (bounds.pMax[dim] - bounds.pMin[dim]);
                for(int i = 0; i < SAH_BUCKETS; i++) {
                    prefix[i] = Bounds3();
                    suffix[i] = Bounds3();
                    buckets[i] = 0;
                }
                
                for(int i = 0; i < objects.size(); i++) {
                    auto cur_bounds = objects[i]->getBounds();
                    double buk_val = (cur_bounds.Centroid()[dim] - bounds.pMin[dim]) * inv_axis_extend * SAH_BUCKETS;
                    int cur_bucket = static_cast<int>(buk_val);
                    assert(cur_bucket >= 0 && cur_bucket < SAH_BUCKETS);
                    buckets[cur_bucket] += 1;
                    prefix[cur_bucket] = Union(prefix[cur_bucket], cur_bounds);
                    suffix[cur_bucket] = Union(suffix[cur_bucket], cur_bounds);
                }

                for(int i = 1; i < SAH_BUCKETS; i++) {
                    prefix[i] = Union(prefix[i], prefix[i - 1]);
                }
                for(int i = SAH_BUCKETS - 2; i >= 0; i--) {
                    suffix[i] = Union(suffix[i], suffix[i + 1]);
                }

                double cur_sah_cost = std::numeric_limits<double>::max();
                int cur_part_bucket = -1;
                for(int i = 0; i < SAH_BUCKETS - 1; i++) {
                    if(buckets[i] == 0) continue;
                    double cur_part_cost = prefix[i].SurfaceArea() * inv_total_area * buckets[i] + suffix[i + 1].SurfaceArea() * inv_total_area * (total_n - buckets[i]) + trav_cost;
                    if(cur_part_cost < cur_sah_cost) {
                        cur_sah_cost = cur_part_cost;
                        cur_part_bucket = i;
                    }
                }
                if(cur_sah_cost < sah_cost) {
                    sah_cost = cur_sah_cost;
                    partion_dim = dim;
                    partion_bucket = cur_part_bucket;
                    assert(cur_part_bucket != -1);
                }
            }

            double leaf_cost = (double)total_n * 1.0;
            if (sah_cost >= leaf_cost || partion_bucket == -1) { 
                node->bounds = bounds;
                node->objects = objects;
                
            } else {
                assert(partion_dim != -1);
                auto leftshapes = std::vector<Object*>();
                auto rightshapes = std::vector<Object*>();
                inv_axis_extend = 1.0 / (bounds.pMax[partion_dim] - bounds.pMin[partion_dim]);

                for(int i = 0; i < objects.size(); i++) {
                    auto cur_bounds = objects[i]->getBounds();
                    double buk_val = (cur_bounds.Centroid()[partion_dim] - bounds.pMin[partion_dim]) * inv_axis_extend * SAH_BUCKETS;
                    int cur_bucket = static_cast<int>(buk_val);
                    assert(partion_bucket != -1);
                    assert(cur_bucket >=0 && cur_bucket < SAH_BUCKETS);
                    if (cur_bucket <= partion_bucket) {
                        leftshapes.push_back(objects[i]);
                    } else {
                        rightshapes.push_back(objects[i]);
                    }
                }
                assert(objects.size() == (leftshapes.size() + rightshapes.size()));

                node->left = recursiveBuild(leftshapes);
                node->right = recursiveBuild(rightshapes);

                node->bounds = Union(node->left->bounds, node->right->bounds);
            }

        }
        break;
    }
    return node;
}

Intersection BVHAccel::Intersect(const Ray& ray) const
{
    Intersection isect;
    if (!root)
        return isect;
    isect = BVHAccel::getIntersection(root, ray);
    return isect;
}

Intersection BVHAccel::getIntersection(BVHBuildNode* node, const Ray& ray) const
{
    // TODO Traverse the BVH to find intersection
    Intersection isect;
    auto dir = ray.direction;
    auto invDir = Vector3f(1.0 / dir.x, 1.0 / dir.y, 1.0 / dir.z);
    std::array<int, 3> dirIsNeg{dir.x > 0, dir.y > 0, dir.z > 0};
    if (node->bounds.IntersectP(ray, invDir, dirIsNeg)) {
        if (!node->left && !node->right) {
            for(auto obj_ptr: node->objects) {
                auto cur_isect = obj_ptr->getIntersection(ray);
                if(cur_isect.distance < isect.distance) isect = cur_isect;
            }
            return isect;
        }
        auto l_isect = BVHAccel::getIntersection(node->left, ray);
        auto r_isect = BVHAccel::getIntersection(node->right, ray);
        isect = (r_isect.distance < l_isect.distance) ? r_isect: l_isect;
    }
    return isect;
}