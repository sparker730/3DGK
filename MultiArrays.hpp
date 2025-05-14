#include <cstddef>

template<typename T>
class Array3D {
public:
    Array3D(){}

    void CreateArray3D(const T* const data, const std::size_t x, const std::size_t y, const std::size_t z) {
       data_ = const_cast<T*>(data);
       x_ = x;
       y_ = y;
       z_ = z;
    }

    int getX(){
        return x_;
    }

    int getY(){
        return y_;
    }

    int getZ(){
        return z_;
    }
    
    inline T& operator()(const std::size_t i, const std::size_t j, const std::size_t k) {
        return data_[k + (j * y_) + (i * x_ * y_)];
    }
private:
    T* data_ = nullptr;
    std::size_t x_ = 0;
    std::size_t y_ = 0;
    std::size_t z_ = 0;
};