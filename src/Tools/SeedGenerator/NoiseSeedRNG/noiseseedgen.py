
def randNoise32(position, seed):
    #Uses bit manipulation to generate 32-bit pseudo-random numbers.

    #Parameters:
    #       position (int): Where in the noise function should be sampled?
    #       seed (int): Input seed for bit mashing.

    #Returns:
    #       mangled (int32): 32-bit bit-mangled value.
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
    #Uses bit manipulation to generate 64-bit pseudo-random numbers.

    #Parameters:
    #       position (int): Where in the noise function should be sampled?
    #       seed (int): Input seed for bit mashing.

    #Returns:
    #       mangled (int64): 64-bit bit-mangled value.

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