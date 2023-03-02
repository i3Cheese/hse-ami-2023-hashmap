/*
 * Cuckoo Hashing
 * 
 * This is a simple implementation of Cuckoo Hashing.
 *
 * Each key has two random buckets in wich it can be stored.
 * On insertion bucket randomly chosen from these two.
 * If the bucket is empty, the key is inserted.
 * If the bucket is not empty, the key is swapped with the key in the bucket.
 * And the process is repeated for the other key.
 * If the number of swaps exceeds a limit, the table is resized, hashes is randomly changes and the process is repeated.
 *
 *
 * Author:  Daniil Sotnikov
 */

#include <algorithm>
#include <cassert>
#include <exception>
#include <initializer_list>
#include <iostream>
#include <limits>
#include <list>
#include <optional>
#include <random>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace __internal_hash_map {

static_assert(sizeof(uint64_t) == 8, "work only for 64bits");
inline const uint64_t __prime_list[] = {
    2ull, 3ull, 5ull, 7ull, 11ull, 13ull, 17ull, 19ull, 23ull, 29ull, 31ull, 37ull, 41ull,
    43ull, 47ull, 53ull, 59ull, 61ull, 67ull, 71ull, 73ull, 79ull, 83ull, 89ull, 97ull,
    103ull, 109ull, 113ull, 127ull, 137ull, 139ull, 149ull, 157ull, 167ull, 179ull, 193ull,
    199ull, 211ull, 227ull, 241ull, 257ull, 277ull, 293ull, 313ull, 337ull, 359ull, 383ull,
    409ull, 439ull, 467ull, 503ull, 541ull, 577ull, 619ull, 661ull, 709ull, 761ull, 823ull,
    887ull, 953ull, 1031ull, 1109ull, 1193ull, 1289ull, 1381ull, 1493ull, 1613ull,
    1741ull, 1879ull, 2029ull, 2179ull, 2357ull, 2549ull, 2753ull, 2971ull, 3209ull,
    3469ull, 3739ull, 4027ull, 4349ull, 4703ull, 5087ull, 5503ull, 5953ull, 6427ull,
    6949ull, 7517ull, 8123ull, 8783ull, 9497ull, 10273ull, 11113ull, 12011ull, 12983ull,
    14033ull, 15173ull, 16411ull, 17749ull, 19183ull, 20753ull, 22447ull, 24281ull,
    26267ull, 28411ull, 30727ull, 33223ull, 35933ull, 38873ull, 42043ull, 45481ull,
    49201ull, 53201ull, 57557ull, 62233ull, 67307ull, 72817ull, 78779ull, 85229ull,
    92203ull, 99733ull, 107897ull, 116731ull, 126271ull, 136607ull, 147793ull,
    159871ull, 172933ull, 187091ull, 202409ull, 218971ull, 236897ull, 256279ull,
    277261ull, 299951ull, 324503ull, 351061ull, 379787ull, 410857ull, 444487ull,
    480881ull, 520241ull, 562841ull, 608903ull, 658753ull, 712697ull, 771049ull,
    834181ull, 902483ull, 976369ull, 1056323ull, 1142821ull, 1236397ull, 1337629ull,
    1447153ull, 1565659ull, 1693859ull, 1832561ull, 1982627ull, 2144977ull, 2320627ull,
    2510653ull, 2716249ull, 2938679ull, 3179303ull, 3439651ull, 3721303ull, 4026031ull,
    4355707ull, 4712381ull, 5098259ull, 5515729ull, 5967347ull, 6456007ull, 6984629ull,
    7556579ull, 8175383ull, 8844859ull, 9569143ull, 10352717ull, 11200489ull,
    12117689ull, 13109983ull, 14183539ull, 15345007ull, 16601593ull, 17961079ull,
    19431899ull, 21023161ull, 22744717ull, 24607243ull, 26622317ull, 28802401ull,
    31160981ull, 33712729ull, 36473443ull, 39460231ull, 42691603ull, 46187573ull,
    49969847ull, 54061849ull, 58488943ull, 63278561ull, 68460391ull, 74066549ull,
    80131819ull, 86693767ull, 93793069ull, 101473717ull, 109783337ull, 118773397ull,
    128499677ull, 139022417ull, 150406843ull, 162723577ull, 176048909ull,
    190465427ull, 206062531ull, 222936881ull, 241193053ull, 260944219ull,
    282312799ull, 305431229ull, 330442829ull, 357502601ull, 386778277ull,
    418451333ull, 452718089ull, 489790921ull, 529899637ull, 573292817ull,
    620239453ull, 671030513ull, 725980837ull, 785430967ull, 849749479ull,
    919334987ull, 994618837ull, 1076067617ull, 1164186217ull, 1259520799ull,
    1362662261ull, 1474249943ull, 1594975441ull, 1725587117ull, 1866894511ull,
    2019773507ull, 2185171673ull, 2364114217ull, 2557710269ull, 2767159799ull,
    2993761039ull, 3238918481ull, 3504151727ull, 3791104843ull, 4101556399ull,
    4294967291ull,
    6442450933ull, 8589934583ull, 12884901857ull, 17179869143ull, 25769803693ull,
    34359738337ull, 51539607367ull, 68719476731ull, 103079215087ull, 137438953447ull,
    206158430123ull, 274877906899ull, 412316860387ull, 549755813881ull,
    824633720731ull, 1099511627689ull, 1649267441579ull, 2199023255531ull,
    3298534883309ull, 4398046511093ull, 6597069766607ull, 8796093022151ull,
    13194139533241ull, 17592186044399ull, 26388279066581ull, 35184372088777ull,
    52776558133177ull, 70368744177643ull, 105553116266399ull, 140737488355213ull,
    211106232532861ull, 281474976710597ull, 562949953421231ull, 1125899906842597ull,
    2251799813685119ull, 4503599627370449ull, 9007199254740881ull,
    18014398509481951ull, 36028797018963913ull, 72057594037927931ull,
    144115188075855859ull, 288230376151711717ull, 576460752303423433ull,
    1152921504606846883ull, 2305843009213693951ull, 4611686018427387847ull,
    9223372036854775783ull, 18446744073709551557ull, 18446744073709551557ull};

inline static const size_t __prime_list_size =
    sizeof(__prime_list) / sizeof(uint64_t);
inline static const uint64_t *__prime_list_end =
    __prime_list + __prime_list_size;

inline static uint64_t get_next_prime(uint64_t n) {
    auto it = std::lower_bound(__prime_list, __prime_list_end, n);
    if (it == __prime_list_end) {
        return *(--it);
    }
    ++it;
    return *it;
}

inline static uint64_t get_random_prime(uint64_t rand_num) {
    return __prime_list[rand_num % __prime_list_size];
}

template <class KeyType, class ValueType, class Hash = std::hash<KeyType>>
class HashMap {
  public:
    using Bucket = std::list<std::pair<const KeyType, ValueType>>;
    using Buckets = std::vector<Bucket>;
    Buckets buckets_;
    size_t size_ = 0;

    double kMaxLoadFactor = 0.20;
    Hash hasher_;
    std::mt19937_64 rand_;

    template <bool _IsConst> class RawIterator {
      private:
        using BucketsIt =
            typename std::conditional<_IsConst,
                                      typename Buckets::const_iterator,
                                      typename Buckets::iterator>::type;
        using BucketIt =
            typename std::conditional<_IsConst, typename Bucket::const_iterator,
                                      typename Bucket::iterator>::type;
        using ReturnType = typename std::conditional<
            _IsConst, const std::pair<const KeyType, ValueType> &,
            std::pair<const KeyType, ValueType> &>::type;

      public:
        RawIterator() : its_(std::nullopt) {}

        RawIterator<_IsConst> &operator++() {
            if (!its_.has_value()) {
                return *this;
            }
            ++(its_->bucket_it);
            while (its_->buckets_it->end() == its_->bucket_it) {
                ++its_->buckets_it;
                if (its_->buckets_it == its_->buckets_end) {
                    its_ = std::nullopt;
                    return *this;
                }
                its_->bucket_it = its_->buckets_it->begin();
            }
            return *this;
        }

        RawIterator<_IsConst> operator++(int) {
            auto old = *this;
            ++(*this);
            return old;
        }

        bool operator==(const RawIterator &other) const {
            return its_ == other.its_;
        };

        bool operator!=(const RawIterator &other) const {
            return !(*this == other);
        };

        ReturnType operator*() { return *(its_->bucket_it); };

        auto operator->() { return its_->bucket_it; }

      private:

        struct Iterators {
            BucketsIt buckets_it;
            BucketsIt buckets_end;
            BucketIt bucket_it;

            bool operator==(const Iterators &other) const {
                return buckets_it == other.buckets_it &&
                       buckets_end == other.buckets_end &&
                       bucket_it == other.bucket_it;
            };

            bool operator!=(const Iterators &other) const {
                return !(*this == other);
            };
        };

        std::optional<Iterators> its_;

        RawIterator(BucketsIt buckets_it, BucketsIt buckets_end,
                    BucketIt bucket_it)
            : its_(Iterators{buckets_it, buckets_end, bucket_it}) {
            if (its_->buckets_it == its_->buckets_end) {
                its_ = std::nullopt;
            } else if (its_->buckets_it->end() == its_->bucket_it) {
                ++(*this);
            }
        }

        friend HashMap;

    };

  public:
    using iterator = RawIterator<false>;
    using const_iterator = RawIterator<true>;

    HashMap(const Hash &hasher = Hash()) : hasher_(hasher) {}

    template <typename It>
    HashMap(It begin, It end, const Hash &hasher = Hash()) : hasher_(hasher) {
        regenerate_coefs();
        for (; begin != end; ++begin) {
            insert(*begin);
        }
    }

    HashMap(std::initializer_list<std::pair<KeyType, ValueType>> list,
            const Hash &hasher = Hash())
        : HashMap(list.begin(), list.end(), hasher) {
    }

    HashMap &operator=(const HashMap &other) {
        if (&other == this) {
            return *this;
        }
        size_ = other.size_;
        rand_num_1_ = other.rand_num_1_;
        rand_num_2_ = other.rand_num_2_;
        hasher_ = other.hasher_;

        // strange copy cause std:pair<const K, V> doesnt have copy assignment operator
        buckets_.resize(other.buckets_.size());
        for (size_t i = 0; i < other.buckets_.size(); ++i) {
            buckets_[i].clear();
            for (const auto &p : other.buckets_[i]) {
                buckets_[i].push_back(p);
            }
        }
        return *this;
    }

    size_t size() const { return size_; }

    bool empty() const { return size_ == 0; }

    Hash hash_function() const { return hasher_; }

    std::pair<iterator, bool>
    insert(const std::pair<KeyType, ValueType> &data) {
        auto it = find(data.first);
        if (it != end()) {
            return {it, false};
        }
        enhance();
        Bucket bucket;
        bucket.push_back(data);
        while (!insert_bucket(bucket)) {
            rehash();
        }
        return {find(data.first), true};
    }

    void erase(const KeyType &key) { erase(find(key)); }

    void erase(const iterator &it) {
        if (it == end()) {
            return;
        }
        it.its_->buckets_it->erase(it.its_->bucket_it);
        --size_;
    }

    iterator begin() {
        if (empty() || buckets_.empty()) {
            return end();
        }
        return iterator(buckets_.begin(), buckets_.end(),
                        buckets_.front().begin());
    }

    iterator end() { return iterator(); }

    const_iterator begin() const {
        if (empty() || buckets_.empty()) {
            return end();
        }
        return const_iterator(buckets_.begin(), buckets_.end(),
                              buckets_.front().begin());
    }

    const_iterator end() const { return const_iterator(); }

    iterator find(const KeyType &key) {
        if (buckets_.empty()) {
            return end();
        }
        {
            const size_t bucket_num = get_hash_1(key) % buckets_.size();
            auto buckets_it = buckets_.begin() + bucket_num;
            for (auto it = buckets_it->begin(); it != buckets_it->end(); ++it) {
                if (it->first == key) {
                    return iterator(buckets_it, buckets_.end(), it);
                }
            }
        }
        { // Differs in GetHash num
            const size_t bucket_num = get_hash_2(key) % buckets_.size();
            auto buckets_it = buckets_.begin() + bucket_num;
            for (auto it = buckets_it->begin(); it != buckets_it->end(); ++it) {
                if (it->first == key) {
                    return iterator(buckets_it, buckets_.end(), it);
                }
            }
        }
        return end();
    }

    const_iterator find(const KeyType &key) const {
        if (buckets_.empty()) {
            return end();
        }
        {
            const size_t bucket_num = get_hash_1(key) % buckets_.size();
            auto buckets_it = buckets_.begin() + bucket_num;
            for (auto it = buckets_it->begin(); it != buckets_it->end(); ++it) {
                if (it->first == key) {
                    return const_iterator(buckets_it, buckets_.end(), it);
                }
            }
        }
        { // Differs in GetHash num
            const size_t bucket_num = get_hash_2(key) % buckets_.size();
            auto buckets_it = buckets_.begin() + bucket_num;
            for (auto it = buckets_it->begin(); it != buckets_it->end(); ++it) {
                if (it->first == key) {
                    return const_iterator(buckets_it, buckets_.end(), it);
                }
            }
        }
        return end();
    }

    ValueType &operator[](const KeyType &key) {
        auto [iter, inserted] = insert({key, ValueType()});
        return iter->second;
    }

    const ValueType &at(const KeyType &key) const {
        auto iter = find(key);
        if (iter == end()) {
            throw std::out_of_range("Not found key");
        }
        std::unordered_map<int, int> a;
        a.begin();
        return iter->second;
    }

    void clear() {
        buckets_.clear();
        size_ = 0;
    }

  private:
    uint64_t rand_num_1_;
    uint64_t rand_num_2_;

    void regenerate_coefs() {
        rand_num_1_ = get_random_prime(rand_());
        rand_num_2_ = get_random_prime(rand_());
    }

    size_t get_hash_1(const KeyType &key) const {
        return hasher_(key) ^ rand_num_1_;
    }

    size_t get_hash_2(const KeyType &key) const {
        return hasher_(key) ^ rand_num_2_;
    }

    void enhance() {
        ++size_;
        if (buckets_.size() == 0) {
            buckets_.resize(get_next_prime(0));
        }
    }

    void rehash() { 
        size_t new_size = std::min((size_t)(kMaxLoadFactor * size_), buckets_.size() * 2);
        new_size = std::max(buckets_.size(), new_size);
        rehash(get_next_prime(new_size));
    }

    void rehash(size_t new_buckets_size) {
        bool success = false;
        while (!success) {
            buckets_.resize(new_buckets_size);
            regenerate_coefs();
            success = true;
            for (size_t i = 0; i < buckets_.size(); ++i) {
                if (!insert_bucket(buckets_[i])) {
                    success = false;
                    break;
                };
            }
            new_buckets_size = get_next_prime(2 * buckets_.size());
        }
    }

    bool insert_bucket(Bucket &to_insert) {
        if (to_insert.empty()) {
            return true;
        }
        const size_t max_cnt_tries = std::sqrt(buckets_.size()) + 10;
        for (size_t cnt_tries = 0; cnt_tries < max_cnt_tries; ++cnt_tries) {
            const KeyType &key = to_insert.front().first;
            size_t bucket_num = (rand_() % 2 == 0)
                                    ? (get_hash_1(key) % buckets_.size())
                                    : (get_hash_2(key) % buckets_.size());
            Bucket &bucket = buckets_[bucket_num];
            if (&bucket == &to_insert) {
                return true;
            }
            // Only hash indistinguishable keys go into same bucket
            if (bucket.empty() ||
                hasher_(bucket.front().first) == hasher_(key)) {
                bucket.splice(bucket.end(), to_insert);
                return true;
            } else {
                to_insert.swap(bucket);
            }
        }
        return false;
    }
};

} // __internal_hash_map

using __internal_hash_map::HashMap;
