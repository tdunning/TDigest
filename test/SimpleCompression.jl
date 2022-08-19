# Licensed to Ted Dunning under one or more
# contributor license agreements.  See the NOTICE file distributed with
# this work for additional information regarding copyright ownership.
# The ASF licenses this file to You under the Apache License, Version 2.0
# (the "License"); you may not use this file except in compliance with
# the License.  You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

# Copyright 2022 Ted Dunning

using DataStructures


abstract type SimpleEncoder end

"""

Very simple variable byte encoding of 60-bit unsigned integers that
always uses 64bit units for output.  The idea is that the several
values are smashed into 64 bits keeping back a few bits to indicate
how the values are fitted in. The rest of the bits hold the compressed
values into equal-sized chunks.

In this encoding, 4 bits out of each 64 byte value are used to
indicate how the remaining 60 bits are divided. The possible values of
the code are shown in the following table:

| Code | Arrangement |
| ---- | ----------- |
|   14 | 1 X 60 bits  |
|   13 | 2 X 30 bits  |
|   12 | 3 X 20 bits  |
|   11 | 4 X 15 bits  |
|   10 | 5 X 12 bits  |
|    9 | 6 X 10 bits  |
|    8 |  7 X 8 bits  |
|    7 |  8 X 7 bits  |
|    6 | 10 X 6 bits  |
|    5 | 12 X 5 bits  |
|    4 | 15 X 4 bits  |
|    3 | 20 X 3 bits  |
|    2 | 30 X 2 bits  |
|    1 | 60 X 1 bits  |
Size codes for Simple64 compression

| Code | Arrangement  |
| ---- | ------------ |
|    1 | 28 X  1 bits |
|    2 | 14 X  2 bits |
|    3 |  9 X  3 bits | 
|    4 |  7 X  4 bits | 
|    5 |  5 X  5 bits | 
|    6 |  4 X  7 bits | 
|    7 |  3 X  9 bits | 
|    8 |  2 X 14 bits |
|    9 |  1 X 28 bits |
Size codes for Simple32 compression

The algorithm for coding a sequence of values is a greedy algorithm
that finds the smallest size code that can encode the next few
values. It does this by starting with the smallest possible size code
for the next value to encode and checks the corresponding range of
upcoming values. If any value in that range is too large for that size
code, the next larger size code is tested with a correspondingly
shorter range of values until a workable code is found.

"""
mutable struct Simple64 <: SimpleEncoder
    values::CircularBuffer
    minCode::CircularBuffer
    maxCode
    function Simple64()
        new(CircularBuffer{UInt64}(61), CircularBuffer{UInt64}(61), 0)
    end
end

mutable struct Simple32 <: SimpleEncoder
    values::CircularBuffer{UInt32}
    minCode::CircularBuffer{UInt32}
    maxCode
    function Simple32()
        new(CircularBuffer{UInt32}(33), CircularBuffer{UInt32}(33), 0)
    end
end

const NUM_DATA_BITS = 60
const BITS_30_MASK = (1 << 30) - 1
const BITS_20_MASK = (1 << 20) - 1
const BITS_15_MASK = (1 << 15) - 1
const BITS_12_MASK = (1 << 12) - 1
const BITS_11_MASK = (1 << 11) - 1
const BITS_10_MASK = (1 << 10) - 1
const BITS_8_MASK = (1 << 8) - 1
const BITS_7_MASK = (1 << 7) - 1
const BITS_6_MASK = (1 << 6) - 1
const BITS_5_MASK = (1 << 5) - 1
const BITS_4_MASK = (1 << 4) - 1
const BITS_3_MASK = (1 << 3) - 1
const BITS_2_MASK = (1 << 2) - 1
const BITS_1_MASK = (1 << 1) - 1

const STATUS_60NUM_1BITS = 1
const STATUS_30NUM_2BITS = 2
const STATUS_20NUM_3BITS = 3
const STATUS_15NUM_4BITS = 4
const STATUS_12NUM_5BITS = 5
const STATUS_10NUM_6BITS = 6
const STATUS_8NUM_7BITS  = 7
const STATUS_7NUM_8BITS  = 8
const STATUS_6NUM_10BITS = 9
const STATUS_5NUM_12BITS = 10
const STATUS_4NUM_15BITS = 11
const STATUS_3NUM_20BITS = 12
const STATUS_2NUM_30BITS = 13
const STATUS_1NUM_60BITS = 14

simple64_sizes = Vector{UInt64}([1, 2, 3, 4, 5, 6, 7, 8, 10, 12, 15, 20, 30, 60])

const STATUS_28NUM_1BITS = 1
const STATUS_14NUM_2BITS = 2
const STATUS_9NUM_3BITS  = 3
const STATUS_7NUM_4BITS  = 4
const STATUS_5NUM_5BITS  = 5
const STATUS_4NUM_7BITS  = 6
const STATUS_3NUM_9BITS  = 7
const STATUS_2NUM_14BITS = 8
const STATUS_1NUM_28BITS = 9

simple32_sizes = Vector{UInt32}([1, 2, 3, 4, 5, 7, 9, 14, 28])

size_map = Dict(Simple32 => simple32_sizes, Simple64 => simple64_sizes)

"""
Returns the minimum encoding required for a value
"""
function bit_code(sizes, v::Integer)
    bits = 8*sizeof(v) - leading_zeros(v)
    i = searchsortedfirst(sizes, bits)
    if i > length(sizes)
        max_bits = maximum(sizes)
        throw(ArgumentError("Value $v has more than $max_bits bits"))
    end
    return i
end

"""
Resets the internal state of an encoder
"""
function reset!(x::SimpleEncoder)
    empty!(x.values)
    empty!(x.minCode)
    x.maxCode = 0
end

"""
    add!(emit::Function, x::Simple64, v::UInt64)

Add a new value to a buffer. If a new output word is generated, call a
function with that output word. Function `emit` may be called more than
once or not all depending on how much data is buffered.

If you want to force all pending output, check `flush!`(@Ref)

# Examples

```julia
s = Simple64()
compressed = Vector{Uint64}()
for i in 1:1000
    add!(s, i) do w
        push!(compressed, w)
    end
end
```

"""
function add!(emit::Function, x::SimpleEncoder, v::Integer)
    sizes = size_map[typeof(x)]
    total_bits = maximum(sizes)

    code = bit_code(sizes, v)
    push!(x.values, v)
    push!(x.minCode, code)
    if code > x.maxCode
        x.maxCode = code
    end

    first_code = x.minCode[1]
    bits = sizes[first_code]
    while sizes[x.maxCode] * length(x.values) ≥ total_bits
        # guaranteed to have at least one output
        n = total_bits ÷ bits
        if n ≤ length(x.values)
            if maximum(x.minCode[1:n]) ≤ first_code
                encode_one(emit, x, first_code)
                break
            end
        end

        first_code += 1
        bits = sizes[first_code]
    end
end

"""
    encode_one(emit, encoder::SimpleEncoder, code)

Emits a full word of compressed data containing values buffered in `encoder`.

Normally only used internally.
"""
function encode_one(emit::Function, x::Simple32, code)
    s9 = code
    if code == STATUS_1NUM_28BITS
        s9 |= popfirst!(x.values) << 4
        popfirst!(x.minCode)

    elseif code == STATUS_2NUM_14BITS
        s9 |= popfirst!(x.values) << 4
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 18
        popfirst!(x.minCode)

    elseif code == STATUS_3NUM_9BITS
        s9 |= popfirst!(x.values) << 4
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 13
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 22
        popfirst!(x.minCode)

    elseif code ==STATUS_4NUM_7BITS
        s9 |= popfirst!(x.values) << 4
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 11
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 18
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 25
        popfirst!(x.minCode)

    elseif code == STATUS_5NUM_5BITS
        s9 |= popfirst!(x.values) << 4
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 9
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 14
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 19
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 24
        popfirst!(x.minCode)

    elseif code ==  STATUS_7NUM_4BITS
        s9 |= popfirst!(x.values) << 4
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 8
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 12
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 16
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 20
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 24
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 28
        popfirst!(x.minCode)

    elseif code ==   STATUS_9NUM_3BITS
        s9 |= popfirst!(x.values) << 4
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 7
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 10
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 13
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 16
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 19
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 22
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 25
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 28
        popfirst!(x.minCode)

    elseif code ==   STATUS_14NUM_2BITS
        s9 |= popfirst!(x.values) << 4
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 6
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 8
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 10
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 12
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 14
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 16
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 18
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 20
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 22
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 24
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 26
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 28
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 30
        popfirst!(x.minCode)

    elseif code ==   STATUS_28NUM_1BITS
        s9 |= popfirst!(x.values)  << 4
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values)  << 5
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values)  << 6
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values)  << 7
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values)  << 8
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values)  << 9
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values)  << 10
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values)  << 11
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values)  << 12
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values)  << 13
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values)  << 14
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values)  << 15
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values)  << 16
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values)  << 17
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values)  << 18
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values)  << 19
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values)  << 20
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values)  << 21
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values)  << 22
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values)  << 23
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values)  << 24
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values)  << 25
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values)  << 26
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values)  << 27
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values)  << 28
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values)  << 29
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values)  << 30
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values)  << 31
        popfirst!(x.minCode)

    else
        throw(ArgumentError("Invalid code: $code"))
    end
    if length(x.minCode) > 0
        x.maxCode = maximum(x.minCode)
    else
        x.maxCode = 0
    end
    emit(s9)
end

function encode_one(emit::Function, x::Simple64, code)
    s9 = code
    if code == STATUS_1NUM_60BITS
        s9 |= popfirst!(x.values) << 4
        popfirst!(x.minCode)

    elseif code == STATUS_2NUM_30BITS
        s9 |= popfirst!(x.values) << 4
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 34
        popfirst!(x.minCode)

    elseif code == STATUS_3NUM_20BITS
        s9 |= popfirst!(x.values) << 4
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 24
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 44
        popfirst!(x.minCode)

    elseif code == STATUS_4NUM_15BITS
        s9 |= popfirst!(x.values) << 4
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 19
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 34
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 49
        popfirst!(x.minCode)

    elseif code == STATUS_5NUM_12BITS
        s9 |= popfirst!(x.values) << 4
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 16
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 28
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 40
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 52
        popfirst!(x.minCode)

    elseif code == STATUS_6NUM_10BITS
        s9 |= popfirst!(x.values) << 4
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 14
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 24
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 34
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 44
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 54
        popfirst!(x.minCode)

    elseif code == STATUS_7NUM_8BITS
        s9 |= popfirst!(x.values) << 4
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 12
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 20
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 28
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 36
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 44
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 52
        popfirst!(x.minCode)

    elseif code == STATUS_8NUM_7BITS
        s9 |= popfirst!(x.values) << 4
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 11
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 18
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 25
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 32
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 39
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 46
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 53
        popfirst!(x.minCode)

    elseif code == STATUS_10NUM_6BITS
        s9 |= popfirst!(x.values) << 4
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 10
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 16
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 22
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 28
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 34
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 40
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 46
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 52
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 58
        popfirst!(x.minCode)

    elseif code == STATUS_12NUM_5BITS
        s9 |= popfirst!(x.values) << 4
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 9
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 14
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 19
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 24
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 29
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 34
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 39
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 44
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 49
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 54
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 59
        popfirst!(x.minCode)

    elseif code == STATUS_15NUM_4BITS
        s9 |= popfirst!(x.values) << 4
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 8
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 12
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 16
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 20
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 24
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 28
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 32
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 36
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 40
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 44
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 48
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 52
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 56
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 60
        popfirst!(x.minCode)

    elseif code == STATUS_20NUM_3BITS
        s9 |= popfirst!(x.values) << 4
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 7
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 10
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 13
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 16
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 19
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 22
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 25
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 28
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 31
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 34
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 37
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 40
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 43
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 46
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 49
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 52
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 55
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 58
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 61
        popfirst!(x.minCode)

    elseif code == STATUS_30NUM_2BITS
        s9 |= popfirst!(x.values) << 4
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 6
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 8
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 10
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 12
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 14
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 16
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 18
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 20
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 22
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 24
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 26
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 28
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 30
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 32
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 34
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 36
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 38
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 40
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 42
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 44
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 46
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 48
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 50
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 52
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 54
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 56
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 58
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 60
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 62
        popfirst!(x.minCode)

    elseif code == STATUS_60NUM_1BITS
        s9 |= popfirst!(x.values) << 4
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 5
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 6
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 7
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 8
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 9
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 10
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 11
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 12
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 13
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 14
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 15
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 16
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 17
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 18
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 19
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 20
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 21
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 22
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 23
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 24
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 25
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 26
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 27
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 28
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 29
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 30
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 31
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 32
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 33
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 34
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 35
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 36
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 37
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 38
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 39
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 40
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 41
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 42
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 43
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 44
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 45
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 46
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 47
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 48
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 49
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 50
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 51
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 52
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 53
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 54
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 55
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 56
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 57
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 58
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 59
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 60
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 61
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 62
        popfirst!(x.minCode)
        s9 |= popfirst!(x.values) << 63
        popfirst!(x.minCode)
    else
        throw(ArgumentError("Invalid code: $code"))
    end
    if length(x.minCode) > 0
        x.maxCode = maximum(x.minCode)
    else
        x.maxCode = 0
    end
    emit(s9)
end

function flush!(emit, x)
    while length(x.values) > 0
        add!(emit, x, 0)
    end
end

function decode_one(emit::Function, input::UInt32)
    code = input & 0xf
    input >>= 4

    if code == STATUS_28NUM_1BITS
        emit((input) & 1)
        emit((input >> 1) & 0x1)
        emit((input >> 2) & 0x1)
        emit((input >> 3) & 0x1)
        emit((input >> 4) & 0x1)
        emit((input >> 5) & 0x1)
        emit((input >> 6) & 0x1)
        emit((input >> 7) & 0x1)
        emit((input >> 8) & 0x1)
        emit((input >> 9) & 0x1)
        emit((input >> 10) & 0x1)
        emit((input >> 11) & 0x1)
        emit((input >> 12) & 0x1)
        emit((input >> 13) & 0x1)
        emit((input >> 14) & 0x1)
        emit((input >> 15) & 0x1)
        emit((input >> 16) & 0x1)
        emit((input >> 17) & 0x1)
        emit((input >> 18) & 0x1)
        emit((input >> 19) & 0x1)
        emit((input >> 20) & 0x1)
        emit((input >> 21) & 0x1)
        emit((input >> 22) & 0x1)
        emit((input >> 23) & 0x1)
        emit((input >> 24) & 0x1)
        emit((input >> 25) & 0x1)
        emit((input >> 26) & 0x1)
        emit((input >> 27) & 0x1)

    elseif code == STATUS_14NUM_2BITS
        emit((input) & 0x3)
        emit((input >> 2) & 0x3)
        emit((input >> 4) & 0x3)
        emit((input >> 6) & 0x3)
        emit((input >> 8) & 0x3)
        emit((input >> 10) & 0x3)
        emit((input >> 12) & 0x3)
        emit((input >> 14) & 0x3)
        emit((input >> 16) & 0x3)
        emit((input >> 18) & 0x3)
        emit((input >> 20) & 0x3)
        emit((input >> 22) & 0x3)
        emit((input >> 24) & 0x3)
        emit((input >> 26) & 0x3)

    elseif code == STATUS_9NUM_3BITS
        emit((input) & 0x7)
        emit((input >> 3) & 0x7)
        emit((input >> 6) & 0x7)
        emit((input >> 9) & 0x7)
        emit((input >> 12) & 0x7)
        emit((input >> 15) & 0x7)
        emit((input >> 18) & 0x7)
        emit((input >> 21) & 0x7)
        emit((input >> 24) & 0x7)

    elseif code == STATUS_7NUM_4BITS
        emit((input) & 0xf)
        emit((input >> 4) & 0xf)
        emit((input >> 8) & 0xf)
        emit((input >> 12) & 0xf)
        emit((input >> 16) & 0xf)
        emit((input >> 20) & 0xf)
        emit((input >> 24) & 0xf)

    elseif code == STATUS_5NUM_5BITS
        emit((input) & 0x1f)
        emit((input >> 5) & 0x1f)
        emit((input >> 10) & 0x1f)
        emit((input >> 15) & 0x1f)
        emit((input >> 20) & 0x1f)

    elseif code == STATUS_4NUM_7BITS
        emit((input) & 0x7f)
        emit((input >> 7) & 0x7f)
        emit((input >> 14) & 0x7f)
        emit((input >> 21) & 0x7f)

    elseif code == STATUS_3NUM_9BITS 
        emit((input) & 0x1ff)
        emit((input >> 9) & 0x1ff)
        emit((input >> 18) & 0x1ff)

    elseif code == STATUS_2NUM_14BITS
        emit((input) & 0x3FFF)
        emit((input >> 14) & 0x3FFF)

    elseif code == STATUS_1NUM_28BITS
        emit(input)

    else
        throw(ArgumentError("Invalid code ($code) in compressed value"))
    end
end

function decode_one(emit::Function, input::UInt64)
    code = input & 0xf
    
    if code == STATUS_1NUM_60BITS
        emit( input >>> 4)

    elseif code ==  STATUS_2NUM_30BITS
        emit((input >>> 4) & BITS_30_MASK)
        emit((input >>> 34) & BITS_30_MASK)

    elseif code ==  STATUS_3NUM_20BITS
        emit((input >>> 4) & BITS_20_MASK)
        emit(((input >>> 24)) & BITS_20_MASK)
        emit((input >>> 44) & BITS_20_MASK)

    elseif code ==  STATUS_4NUM_15BITS
        emit((input >>> 4) & BITS_15_MASK)
        emit((input >>> 19) & BITS_15_MASK)
        emit((input >>> 34) & BITS_15_MASK)
        emit((input >>> 49) & BITS_15_MASK)

    elseif code ==  STATUS_5NUM_12BITS
        emit((input >>> 4) & BITS_12_MASK)
        emit((input >>> 16) & BITS_12_MASK)
        emit((input >>> 28) & BITS_12_MASK)
        emit((input >>> 40) & BITS_12_MASK)
        emit((input >>> 52) & BITS_12_MASK)

    elseif code ==  STATUS_6NUM_10BITS
        emit((input >>> 4) & BITS_10_MASK)
        emit((input >>> 14) & BITS_10_MASK)
        emit((input >>> 24) & BITS_10_MASK)
        emit((input >>> 34) & BITS_10_MASK)
        emit((input >>> 44) & BITS_10_MASK)
        emit((input >>> 54) & BITS_10_MASK)

    elseif code ==  STATUS_7NUM_8BITS
        emit((input >>> 4) & BITS_8_MASK)
        emit((input >>> 12) & BITS_8_MASK)
        emit((input >>> 20) & BITS_8_MASK)
        emit((input >>> 28) & BITS_8_MASK)
        emit((input >>> 36) & BITS_8_MASK)
        emit((input >>> 44) & BITS_8_MASK)
        emit((input >>> 52) & BITS_12_MASK)

    elseif code ==  STATUS_8NUM_7BITS
        emit((input >>> 4) & BITS_7_MASK)
        emit((input >>> 11) & BITS_7_MASK)
        emit((input >>> 18) & BITS_7_MASK)
        emit((input >>> 25) & BITS_7_MASK)
        emit((input >>> 32) & BITS_7_MASK)
        emit((input >>> 39) & BITS_7_MASK)
        emit((input >>> 46) & BITS_7_MASK)
        emit((input >>> 53) & BITS_11_MASK)

    elseif code ==  STATUS_10NUM_6BITS
        emit((input >>> 4) & BITS_6_MASK)
        emit((input >>> 10) & BITS_6_MASK)
        emit((input >>> 16) & BITS_6_MASK)
        emit((input >>> 22) & BITS_6_MASK)
        emit((input >>> 28) & BITS_6_MASK)
        emit((input >>> 34) & BITS_6_MASK)
        emit((input >>> 40) & BITS_6_MASK)
        emit((input >>> 46) & BITS_6_MASK)
        emit((input >>> 52) & BITS_6_MASK)
        emit((input >>> 58) & BITS_6_MASK)

    elseif code ==  STATUS_12NUM_5BITS
        emit((input >>> 4) & BITS_5_MASK)
        emit((input >>> 9) & BITS_5_MASK)
        emit((input >>> 14) & BITS_5_MASK)
        emit((input >>> 19) & BITS_5_MASK)
        emit((input >>> 24) & BITS_5_MASK)
        emit((input >>> 29) & BITS_5_MASK)
        emit((input >>> 34) & BITS_5_MASK)
        emit((input >>> 39) & BITS_5_MASK)
        emit((input >>> 44) & BITS_5_MASK)
        emit((input >>> 49) & BITS_5_MASK)
        emit((input >>> 54) & BITS_5_MASK)
        emit((input >>> 59) & BITS_5_MASK)

    elseif code ==  STATUS_15NUM_4BITS
        emit((input >>> 4) & BITS_4_MASK)
        emit((input >>> 8) & BITS_4_MASK)
        emit((input >>> 12) & BITS_4_MASK)
        emit((input >>> 16) & BITS_4_MASK)
        emit((input >>> 20) & BITS_4_MASK)
        emit((input >>> 24) & BITS_4_MASK)
        emit((input >>> 28) & BITS_4_MASK)
        emit((input >>> 32) & BITS_4_MASK)
        emit((input >>> 36) & BITS_4_MASK)
        emit((input >>> 40) & BITS_4_MASK)
        emit((input >>> 44) & BITS_4_MASK)
        emit((input >>> 48) & BITS_4_MASK)
        emit((input >>> 52) & BITS_4_MASK)
        emit((input >>> 56) & BITS_4_MASK)
        emit((input >>> 60) & BITS_4_MASK)

    elseif code ==  STATUS_20NUM_3BITS
        emit((input >>> 4) & BITS_3_MASK)
        emit((input >>> 7) & BITS_3_MASK)
        emit((input >>> 10) & BITS_3_MASK)
        emit((input >>> 13) & BITS_3_MASK)
        emit((input >>> 16) & BITS_3_MASK)
        emit((input >>> 19) & BITS_3_MASK)
        emit((input >>> 22) & BITS_3_MASK)
        emit((input >>> 25) & BITS_3_MASK)
        emit((input >>> 28) & BITS_3_MASK)
        emit((input >>> 31) & BITS_3_MASK)
        emit((input >>> 34) & BITS_3_MASK)
        emit((input >>> 37) & BITS_3_MASK)
        emit((input >>> 40) & BITS_3_MASK)
        emit((input >>> 43) & BITS_3_MASK)
        emit((input >>> 46) & BITS_3_MASK)
        emit((input >>> 49) & BITS_3_MASK)
        emit((input >>> 52) & BITS_3_MASK)
        emit((input >>> 55) & BITS_3_MASK)
        emit((input >>> 58) & BITS_3_MASK)
        emit((input >>> 61) & BITS_3_MASK)

    elseif code ==  STATUS_30NUM_2BITS
        emit((input >>> 4) & BITS_2_MASK)
        emit((input >>> 6) & BITS_2_MASK)
        emit((input >>> 8) & BITS_2_MASK)
        emit((input >>> 10) & BITS_2_MASK)
        emit((input >>> 12) & BITS_2_MASK)
        emit((input >>> 14) & BITS_2_MASK)
        emit((input >>> 16) & BITS_2_MASK)
        emit((input >>> 18) & BITS_2_MASK)
        emit((input >>> 20) & BITS_2_MASK)
        emit((input >>> 22) & BITS_2_MASK)
        emit((input >>> 24) & BITS_2_MASK)
        emit((input >>> 26) & BITS_2_MASK)
        emit((input >>> 28) & BITS_2_MASK)
        emit((input >>> 30) & BITS_2_MASK)
        emit((input >>> 32) & BITS_2_MASK)
        emit((input >>> 34) & BITS_2_MASK)
        emit((input >>> 36) & BITS_2_MASK)
        emit((input >>> 38) & BITS_2_MASK)
        emit((input >>> 40) & BITS_2_MASK)
        emit((input >>> 42) & BITS_2_MASK)
        emit((input >>> 44) & BITS_2_MASK)
        emit((input >>> 46) & BITS_2_MASK)
        emit((input >>> 48) & BITS_2_MASK)
        emit((input >>> 50) & BITS_2_MASK)
        emit((input >>> 52) & BITS_2_MASK)
        emit((input >>> 54) & BITS_2_MASK)
        emit((input >>> 56) & BITS_2_MASK)
        emit((input >>> 58) & BITS_2_MASK)
        emit((input >>> 60) & BITS_2_MASK)
        emit((input >>> 62) & BITS_2_MASK)

    elseif code ==  STATUS_60NUM_1BITS
        emit((input >>> 4) & BITS_1_MASK)
        emit((input >>> 5) & BITS_1_MASK)
        emit((input >>> 6) & BITS_1_MASK)
        emit((input >>> 7) & BITS_1_MASK)
        emit((input >>> 8) & BITS_1_MASK)
        emit((input >>> 9) & BITS_1_MASK)
        emit((input >>> 10) & BITS_1_MASK)
        emit((input >>> 11) & BITS_1_MASK)
        emit((input >>> 12) & BITS_1_MASK)
        emit((input >>> 13) & BITS_1_MASK)
        emit((input >>> 14) & BITS_1_MASK)
        emit((input >>> 15) & BITS_1_MASK)
        emit((input >>> 16) & BITS_1_MASK)
        emit((input >>> 17) & BITS_1_MASK)
        emit((input >>> 18) & BITS_1_MASK)
        emit((input >>> 19) & BITS_1_MASK)
        emit((input >>> 20) & BITS_1_MASK)
        emit((input >>> 21) & BITS_1_MASK)
        emit((input >>> 22) & BITS_1_MASK)
        emit((input >>> 23) & BITS_1_MASK)
        emit((input >>> 24) & BITS_1_MASK)
        emit((input >>> 25) & BITS_1_MASK)
        emit((input >>> 26) & BITS_1_MASK)
        emit((input >>> 27) & BITS_1_MASK)
        emit((input >>> 28) & BITS_1_MASK)
        emit((input >>> 29) & BITS_1_MASK)
        emit((input >>> 30) & BITS_1_MASK)
        emit((input >>> 31) & BITS_1_MASK)
        emit((input >>> 32) & BITS_1_MASK)
        emit((input >>> 33) & BITS_1_MASK)
        emit((input >>> 34) & BITS_1_MASK)
        emit((input >>> 35) & BITS_1_MASK)
        emit((input >>> 36) & BITS_1_MASK)
        emit((input >>> 37) & BITS_1_MASK)
        emit((input >>> 38) & BITS_1_MASK)
        emit((input >>> 39) & BITS_1_MASK)
        emit((input >>> 40) & BITS_1_MASK)
        emit((input >>> 41) & BITS_1_MASK)
        emit((input >>> 42) & BITS_1_MASK)
        emit((input >>> 43) & BITS_1_MASK)
        emit((input >>> 44) & BITS_1_MASK)
        emit((input >>> 45) & BITS_1_MASK)
        emit((input >>> 46) & BITS_1_MASK)
        emit((input >>> 47) & BITS_1_MASK)
        emit((input >>> 48) & BITS_1_MASK)
        emit((input >>> 49) & BITS_1_MASK)
        emit((input >>> 50) & BITS_1_MASK)
        emit((input >>> 51) & BITS_1_MASK)
        emit((input >>> 52) & BITS_1_MASK)
        emit((input >>> 53) & BITS_1_MASK)
        emit((input >>> 54) & BITS_1_MASK)
        emit((input >>> 55) & BITS_1_MASK)
        emit((input >>> 56) & BITS_1_MASK)
        emit((input >>> 57) & BITS_1_MASK)
        emit((input >>> 58) & BITS_1_MASK)
        emit((input >>> 59) & BITS_1_MASK)
        emit((input >>> 60) & BITS_1_MASK)
        emit((input >>> 61) & BITS_1_MASK)
        emit((input >>> 62) & BITS_1_MASK)
        emit((input >>> 63) & BITS_1_MASK)

    else
        throw(ArgumentError("Unknown code: $code"))
    end
end

"""
    compress!(compressed, uncompressed)

Compresses the values in `uncompressed` into `compressed`. The values
may be padded with zeros so a round trip compression/uncompression of
`n` values will agree in the first `n` values, but may have appended
zeros.

"""
function compress!(compressed::Vector{UInt64}, uncompressed::Vector{<:Integer})
    simple = Simple64()
    for x in uncompressed
        add!(simple, x) do z
            push!(compressed, z)
        end
    end
    flush!(simple) do z
        push!(compressed, z)
    end
    return compressed
end

function compress!(compressed::Vector{UInt32}, uncompressed::Vector{<:Integer})
    simple = Simple32()
    for x in uncompressed
        add!(simple, x) do z
            push!(compressed, z)
        end
    end
    flush!(simple) do z
        push!(compressed, z)
    end
    return compressed
end

function uncompress!(uncompressed, compressed::Vector)
    for x in compressed
        decode_one(x) do v
            push!(uncompressed, v)
        end
    end
    return uncompressed
end

function uncompress!(uncompressed, compressed::Vector, n)
    for x in @view(compressed[begin:(begin+n-1)])
        uncompress(uncompressed, x)
    end
    return uncompressed
end
