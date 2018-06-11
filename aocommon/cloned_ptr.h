#ifndef CLONED_PTR_H
#define CLONED_PTR_H

#include <memory>

namespace ao {

template<typename T>
class cloned_ptr {
public:
	cloned_ptr() noexcept {}
	cloned_ptr(nullptr_t) noexcept {}
	cloned_ptr(T* object) noexcept : _ptr(object) {}
	cloned_ptr(const cloned_ptr<T>& other) :
		_ptr(other._ptr==nullptr ? nullptr : new T(*other._ptr))
	{	}
	cloned_ptr(cloned_ptr<T>&& other) noexcept :
		_ptr(std::move(other._ptr))
	{ }
	
	cloned_ptr<T>& operator=(nullptr_t) noexcept
	{
		_ptr.reset();
	}
	cloned_ptr<T>& operator=(const cloned_ptr<T>& other)
	{
		_ptr.reset(other._ptr==nullptr ? nullptr : new T(*other._ptr));
	}
	cloned_ptr<T>& operator=(cloned_ptr<T>&& other) noexcept
	{
		_ptr = std::move(other._ptr);
	}
	void reset(T* object) noexcept  
	{
		_ptr.reset(object);
	}
	
	T& operator*() const noexcept { return *_ptr; }
	T* operator->() const noexcept { return _ptr.get(); }
	T* get() const { return _ptr.get(); }
	
	bool operator==(nullptr_t) const noexcept
	{
		return _ptr == nullptr;
	}
	bool operator==(const std::unique_ptr<T>& rhs) const noexcept
	{
		return _ptr == rhs;
	}
	bool operator==(const cloned_ptr<T>& rhs) const noexcept
	{
		return _ptr == rhs._ptr;
	}
	void swap(cloned_ptr<T>& other) noexcept
	{
		std::swap(_ptr, other._ptr);
	}
private:
	std::unique_ptr<T> _ptr;
};

template<typename T>
void swap(cloned_ptr<T>& left, cloned_ptr<T>& right)
{
	left.swap(right);
}

}

#endif
