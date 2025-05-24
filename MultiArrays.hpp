#include <cstddef>

template<typename T>
class Array3D {
public:
    Array3D(size_t xsize, size_t ysize, size_t zsize) : x_(xsize), y_(ysize), z_(zsize){
        this->data_ = new T[x_ * y_ * z_];
        std::fill(this->data_, this->data_ + (x_*y_*z_), 0);
        size_ = x_*y_*z_;
    }

    ~Array3D() {
        this->x_=0;
        this->y_=0;
        this->z_=0;
        delete[] this->data_;
    }

    void CreateArray3D(const T* const data, const std::size_t x, const std::size_t y, const std::size_t z) { //legacy?
       data_ = const_cast<T*>(data);
       x_ = x;
       y_ = y;
       z_ = z;
       size_ = x_ * y_ * z_;
    }

    void clear() {
        std::fill(data_, data_+size_, T{});//T{} instead of 0 because of complex types where 0 is real
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

    T* start() {
        return data_;
    }

    inline T& operator()(const std::size_t i, const std::size_t j, const std::size_t k) {
        return data_[k + (j * y_) + (i * x_ * y_)];
    }

    inline Array3D& operator=(const Array3D &arr){
        for(auto i = 0; i < size_; ++i){
            data_[i] = arr.data_[i];
        }
        return *this;
    }


private:
    T* data_ = nullptr;
    std::size_t size_ = 0;
    std::size_t x_ = 0;
    std::size_t y_ = 0;
    std::size_t z_ = 0;
};

template<typename T>
class Array2D{
    public:
    Array2D(size_t xsize, size_t ysize) : x_(xsize), y_(ysize){
        this->data_ = new T[x_ * y_];
        std::fill(this->data_, this->data_ + (x_*y_));
        size_ = x_*y_;
    }

    ~Array2D(){
        x_ = 0;
        y_ = 0;
        delete[] this->data_;
    }

    void clear() {
        std::fill(data_, data_+size_, T{});
    }

    T* start() {
        return data_;
    }

    inline T& operator()(const std::size_t i, const std::size_t j) {
        return data_[j * x_ + i];
    }

    inline Array2D& operator=(const Array2D &arr){
        for(auto i = 0; i < size_; ++i){
            data_[i] = arr.data_[i];
        }
        return *this;
    }

    private:
    T* data_ = nullptr;
    size_t x_ = 0;
    size_t y_ = 0;
    size_t size_ = 0;
};
