#ifndef OLD_HADOOP_HPP
#define OLD_HADOOP_HPP

#include <cassert>
#include <cstdint>
#include <fstream>
#include <string>
#include <vector>

namespace lemon {
namespace deprecated {

class Hadoop {
   public:
    Hadoop(std::istream& stream) : stream_(stream) { initialize_(); }

    bool has_next() { return stream_.peek() != std::char_traits<char>::eof(); }

    std::pair<std::vector<char>, std::vector<char>> next() { return read(); }

   private:
    std::istream& stream_;
    std::string marker_ = "";
    std::vector<char> key_;

    /* Uncomment if a more generic solution is desired....
    Loosly based on
    https://github.com/azat-archive/hadoop-io-sequence-reader/blob/master/src/reader.cpp,
    but entirely rewritten
    static uint32_t decode(int8_t ch) {
        if (ch >= -112) {
            return 1;
        } else if (ch < -120) {
            return -119 - ch;
        }
        return -111 - ch;
    }

    static int64_t read(const char *pos, uint32_t &len) {
        if (*pos >= (char)-112) {
            len = 1;
            return *pos;
        } else {
            return read_helper(pos, len);
        }
    }

    static int64_t read_helper(const char *pos, uint32_t &len) {
        bool neg = *pos < -120;
        len = neg ? (-119 - *pos) : (-111 - *pos);
        const char *end = pos + len;
        int64_t value = 0;
        while (++pos < end) {
            value = (value << 8) | *(uint8_t *)pos;
        }
        return neg ? (value ^ -1LL) : value;
    }

    int64_t read_long() {
        char buff[10];
        stream_.read(buff, 1);
        auto len = decode(*buff);
        if (len > 1) {
            !stream_.read(buff + 1, len - 1);
        }
        return read(buff, len);
    }

    void read_string(char *buffer) {
        auto len = read_long();
        stream_.read(buffer, len);
    }

    void initialize() {
        stream_.exceptions(std::ifstream::badbit | std::ifstream::failbit);

        char buffer[1024];

        // Reads the header
        stream_.read(buffer, 4);

        // Read the Key class
        read_string(buffer);

        // Read the Value class
        read_string(buffer);

        auto valueCompression = stream_.get();
        auto blockCompression = stream_.get();

        if (valueCompression != 0 || blockCompression != 0) {
            throw std::runtime_error("Compression not supported");
        }

        int pairs = read_int();

        if (pairs < 0 || pairs > 1024) {
            throw std::runtime_error("Invalid pair count");
        }

        for (size_t i = 0; i < pairs; ++i) {
            // Ignore the metadata
            read_string(buffer);
            read_string(buffer);
        }

        // Read marker
        stream_.read(buffer, 16);

        marker_ = std::string(buffer, 16);
    } */

    void initialize_() {
        // Completly skip the header as it is the same in all RCSB Hadoop files
        stream_.exceptions(std::ifstream::badbit | std::ifstream::failbit);
        char buffer[86];
        stream_.read(buffer, 87);
    }

    int read_int() {
        int ret;
        stream_.read(reinterpret_cast<char*>(&ret), 4);
        return ntohl(ret);
    }

    std::pair<std::vector<char>, std::vector<char>> read() {
        auto sync_check = read_int();

        if (sync_check == -1) {
            std::vector<char> marker(16);
            stream_.read(marker.data(), 16);
            // Only valid if using the full version
            // assert(std::string(marker.data(), 16) == marker_);
            return this->read();
        }

        auto key_length = read_int();

        // Do not check this during runtime as it should all be the same
        assert(key_length >= 4);

        std::vector<char> key(key_length);
        stream_.read(key.data(), key_length);

        assert(sync_check >= 8);

        // Remove junk characters added by Java Serialization
        char junk[4];
        stream_.read(junk, 4);

        int value_length = sync_check - key_length;
        std::vector<char> value(value_length - 4);
        stream_.read(value.data(), value_length - 4);

        return {key, value};
    }
};
}
}

#endif
