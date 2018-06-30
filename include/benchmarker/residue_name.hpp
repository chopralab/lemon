#ifndef RESIDUE_NAME_HPP
#define RESIDUE_NAME_HPP

#include <array>
#include <string>
#include <unordered_map>

namespace benchmarker {

    class ResidueName : private std::array<char, 3> {
        using super = std::array<char, 3>;
        friend bool operator==(const ResidueName& lhs, const ResidueName& rhs);
        friend bool operator==(const ResidueName& lhs, const std::string& rhs);
        friend bool operator==(const std::string& lhs, const ResidueName& rhs);
        friend bool operator!=(const ResidueName& lhs, const ResidueName& rhs);
        friend bool operator!=(const ResidueName& lhs, const std::string& rhs);
        friend bool operator!=(const std::string& lhs, const ResidueName& rhs);
        friend std::ostream& operator<<(std::ostream& os, const ResidueName&);
        using super::operator[];
    public:
        ResidueName(const std::string& s);
        ResidueName(const char* s);

        unsigned short hash() const;
        const std::array<char, 3>& operator*() const;
    };

    struct ResidueNameHash {
        unsigned short operator()(const ResidueName& resn) const;
    };

    inline bool operator==(const ResidueName& lhs, const ResidueName& rhs) {
        return
            (lhs[0] == rhs[0]) &&
            (lhs[1] == rhs[1]) &&
            (lhs[2] == rhs[2]);
    }

    inline bool operator==(const ResidueName& lhs, const std::string& rhs) {
        switch (rhs.length()) {
            case 3:
                return (lhs[0] == rhs[0]) && (lhs[1] == rhs[1]) && (lhs[2] == rhs[2]);
                break;
            case 2:
                return (lhs[0] == rhs[0]) && (lhs[1] == rhs[1]) && (lhs[2] == 0);
                break;
            case 1:
                return (lhs[0] == rhs[0]) && (lhs[1] == 0) && (lhs[2] == 0);
                break;
            default:
                return false;
                break;
        }

        return false;
    }

    inline bool operator==(const std::string& lhs, const ResidueName& rhs) {
        return rhs == lhs;
    }

    inline bool operator!=(const ResidueName& lhs, const ResidueName& rhs) {
        return !(lhs == rhs);
    }

    inline bool operator!=(const ResidueName& lhs, const std::string& rhs) {
        return !(lhs == rhs);
    }

    inline bool operator!=(const std::string& lhs, const ResidueName& rhs) {
        return !(lhs == rhs);
    }

    typedef std::unordered_map<ResidueName, std::size_t, ResidueNameHash> ResidueNameCount;

    ResidueNameCount& operator+=(ResidueNameCount& lhs, const ResidueNameCount& rhs);
    std::ostream& operator<<(std::ostream& os, const ResidueNameCount& rnc);
}

#endif
