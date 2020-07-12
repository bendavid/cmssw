#ifndef DataFormats_Common_RefItemGet_h
#define DataFormats_Common_RefItemGet_h

/*----------------------------------------------------------------------
  
RefItemGet: Free function to get pointer to a referenced item.


----------------------------------------------------------------------*/
#include "DataFormats/Common/interface/RefCore.h"
#include "DataFormats/Common/interface/RefCoreGet.h"

namespace edm {

  namespace refitem {

    template <typename C, typename T, typename F, typename K>
    inline void findRefItem(RefCore const& refCore, C const* container, K const& key) {
      F finder;
      T const* item = finder(*container, key);
      refCore.setProductPtr(item);
    }

    template <typename C, typename T, typename F, typename KEY>
    struct GetRefPtrImpl {
      static T const* getRefPtr_(RefCore const& product, KEY const& key) {
        T const* item = static_cast<T const*>(product.productPtr());
        if (item != nullptr) {
          return item;
        }
        auto prodGetter = product.productGetter();
        if (nullptr == prodGetter) {
          item = static_cast<T const*>(product.productPtr());
          if (item != nullptr) {
            //Another thread updated the value since we checked
            return item;
          }
        }
        C const* prod = edm::template getProductWithCoreFromRef<C>(product, prodGetter);
        /*
        typename C::const_iterator it = prod->begin();
         std::advance(it, item.key());
         T const* p = it.operator->();
        */
        F func;
        item = func(*prod, key);
        product.setProductPtr(item);
        return item;
      }
    };

    template <typename C, typename T, typename F>
    struct GetRefPtrImpl<C, T, F, unsigned int> {
      static T const* getRefPtr_(RefCore const& product, unsigned int key) {
        T const* item = static_cast<T const*>(product.productPtr());
        if (item != nullptr) {
          return item;
        }
        auto getter = product.productGetter();
        if (getter == nullptr) {
          auto prod = product.productPtr();
          if (prod != nullptr) {
            //another thread updated the value since we last checked.
            return static_cast<T const*>(prod);
          }
        }
        C const* prod = edm::template tryToGetProductWithCoreFromRef<C>(product, getter);
        if (prod != nullptr) {
          F func;
          item = func(*prod, key);
          product.setProductPtr(item);
          return item;
        }
        unsigned int thinnedKey = key;
        prod = edm::template getThinnedProduct<C>(product, thinnedKey, getter);
        F func;
        item = func(*prod, thinnedKey);
        product.setProductPtr(item);
        return item;
      }
    };
  }  // namespace refitem

  template <typename C, typename T, typename F, typename KEY>
  inline T const* getRefPtr(RefCore const& product, KEY const& iKey) {
    return refitem::GetRefPtrImpl<C, T, F, KEY>::getRefPtr_(product, iKey);
  }

  namespace refitem {
    template <typename C, typename KEY>
    struct IsThinnedAvailableImpl {
      static bool isThinnedAvailable_(RefCore const& product, KEY const& key) { return false; }
    };

    template <typename C>
    struct IsThinnedAvailableImpl<C, unsigned int> {
      static bool isThinnedAvailable_(RefCore const& ref, unsigned int key) {
        if (ref.productPtr() != nullptr) {
          return true;
        }
        if (ref.isTransient()) {
          return false;
        }
        auto getter = ref.productGetter();
        if (getter != nullptr) {
          return ref.isThinnedAvailable(key, getter);
        }
        //another thread may have updated the cache
        return nullptr != ref.productPtr();
      }
    };
  }  // namespace refitem

  template <typename C, typename KEY>
  inline bool isThinnedAvailable(RefCore const& product, KEY const& iKey) {
    return refitem::IsThinnedAvailableImpl<C, KEY>::isThinnedAvailable_(product, iKey);
  }

  template <typename C>
  inline Ref<C> thinnedRefFrom(Ref<C> const& parent, ProductID const& targetpid, EDProductGetter const& prodGetter) {
    printf("thinnedRefFrom\n");
    if (parent.id() == targetpid) {
      printf("returning parent\n");
      return parent;
    }
    printf("checking for thinned collections\n");
    unsigned int key = parent.key();
    WrapperBase const* thinnedprod = prodGetter.getThinnedProduct(parent.id(), key, targetpid);
    printf("checking result\n");
    if (thinnedprod == nullptr) {
      printf("returning null\n");
      return Ref<C>(targetpid);
    }
    printf("getting item pointer\n");
    typename edm::Ref<C>::finder_type finder;
    Wrapper<C> const* wrapper = static_cast<Wrapper<C> const*>(thinnedprod);
    typename edm::Ref<C>::value_type const* item = finder(*wrapper->product(), key);
    //     typename edm::Ref<C>::value_type const* item =
    //             finder(*(static_cast<typename edm::Ref<C>::product_type const*>(thinnedprod->product())), key);
    printf("returning final ref\n");
    return Ref<C>(targetpid, item, key);
    //TODO we already have the collection pointer, so we could fill the cache of the Ref directly
    //to speed up any subsequent access
    //     return Ref<C>(targetpid, key, &prodGetter);
  }

}  // namespace edm

#endif
