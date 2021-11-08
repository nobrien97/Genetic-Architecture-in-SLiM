
# Generate some random numbers
def randNoise32(position, seed):
    PRIME1 = 0x75BD0FB
    PRIME2 = 0x75BD12D
    mangled = position
    mangled *= PRIME1
    mangled += seed
    mangled ^= (mangled >> 8)
    mangled *= PRIME2
    mangled ^= (mangled << 6)
    mangled *= mangled
    mangled ^= (mangled >> 5)
    # Cast into 32 bit value
    return mangled & 0xffffffff


def randNoise64(position, seed):
    PRIME1 = 0x75BD0FB
    PRIME2 = 0x75BD12D
    mangled = position
    mangled *= PRIME1
    mangled += seed
    mangled ^= (mangled >> 8)
    mangled *= PRIME2
    mangled ^= (mangled << 6)
    mangled *= mangled
    mangled ^= (mangled >> 5)
    # Cast into 64 bit value
    return mangled & 0xffffffffffffffff