#ifndef _LRU_CACHE_H_
#define _LRU_CACHE_H_

#include <list>
#include <map>
#include <optional>
#include <utility>


/**
 * This is an owning least recently used cache. When an item is inserted
 * into the cache, the cache takes ownership of the object. When the
 * cache limit is reached, the least recently used object is evicted from
 * the cache. This will destroy the object.
 * 
 * Because the cache will destroy objects that are evicted, this is best
 * for objects that are expensive to create, but also will take up a lot
 * of memory if the entire collection of objects are held in memory.
 */
template <class Key, class Value>
class LruCache
{
public:

  using KeyType = Key;
  using ValueType = Value;
  using ListType = std::list<KeyType>;
  using MapValueType = std::pair<ValueType, typename ListType::iterator>;
  using MapType = std::map<KeyType, MapValueType>;
  using OptionType = std::optional<ValueType&>;

  LruCache(size_t capacity)
    : m_capacity(capacity)
  {
  }

  // cannot copy
  LruCache(LruCache const& o) = delete;
  LruCache& operator=(LruCache const& o) = delete;

  // moving is allowed, but it's complexity is 
  // on the order of the number of elements in
  // the cache.
  LruCache(LruCache && o);
  LruCache& operator=(LruCache && o);

  ~LruCache() {}

  size_t size() const { return m_map.size(); }

  size_t capacity() const { return m_capacity; }

  bool empty() const { return m_map.empty(); }

  bool contains(const KeyType &key) { return m_map.find(key) != m_map.end(); }

  void insert(const KeyType &key, ValueType value)
  {
    auto i = m_map.find(key);
    if (i == m_map.end())
    {
      // insert item into the cache, but first check if it is full
      if (size() >= m_capacity)
      {
        // cache is full, evict the least recently used item
        evict();
      }

      // insert the new item
      m_list.push_front(key);
      m_map.insert({key, MapValueType(std::move(value), m_list.begin())});
    }
  }

  ValueType & get(KeyType const &key)
  {
    // lookup value in the cache
    auto i = m_map.find(key);
    if (i == m_map.end())
    {
      throw std::runtime_error("No value found for key");
    }

    // return the value, but first update its place in the most
    // recently used list
    auto j = i->second.second;
    if (j != m_list.begin())
    {
      // move item to the front of the most recently used list
      m_list.erase(j);
      m_list.push_front(key);

      // update iterator in map
      j = m_list.begin();
      ValueType &value = i->second.first;
      i->second.second = j;
      // m_map[key] = std::make_pair(value, j);

      // return the value
      return value;
    }
    else
    {
      // the item is already at the front of the most recently
      // used list so just return it
      return i->second.first;
    }
  }

  void clear()
  {
    m_map.clear();
    m_list.clear();
  }

private:
  void evict()
  {
    // evict item from the end of most recently used list
    auto i = --m_list.end();
    m_map.erase(*i);
    m_list.erase(i);
  }

private:
  MapType m_map;
  ListType m_list;
  size_t m_capacity;
};


template <class Key, class Value>
LruCache<Key, Value>::LruCache(LruCache && o)
{
  *this = std::move(o);
}


template <class Key, class Value>
LruCache<Key, Value>&
LruCache<Key, Value>::operator=(LruCache && o)
{
  m_capacity = o.m_capacity;

  // have to loop through the map
  // and add each element to get the
  // ordering correct
  for (auto &kv: o.m_map)
    insert(kv.first, std::move(kv.second.first));
  
  return *this;
}

#endif // _LRU_CACHE_H_