###############################################################################
#                             Types compatibilites                            #
###############################################################################

## Add float arithmetic compatibility with rational type
Base.Float64(x::Nemo.fmpq) = Base.Float64(Rational(x))
Base.Float32(x::Nemo.fmpq) = Base.Float32(Rational(x))
Base.Float16(x::Nemo.fmpq) = Base.Float16(Rational(x))
Base.BigFloat(x::Nemo.fmpq) = Base.BigFloat(Rational(x))

## In case we want to push! a fmpq coeff to a float array
Base.convert(::Type{T}, x::Nemo.fmpq) where {T <: AbstractFloat} = T(Rational(x))
Base.convert(::Type{T}, x::Nemo.fmpz) where {T <: AbstractFloat} = T(Int(x))

## From float to fmpq
(a::Nemo.FlintRationalField)(x::T) where {T <: AbstractFloat} = Nemo.fmpq(rationalize(x))
Nemo.fmpq(x::T) where {T <: AbstractFloat} = Nemo.fmpq(rationalize(x))

## From Nemo integers to Juila rationals
Rational{T}(x::Nemo.fmpz) where {T <: Integer} = Rational{T}(T(x))

## From Polymake matrices to Juila arrays
function convert_to_array(M::Polymake.Matrix, ::T=1) where {T}
    mat = Polymake.@convert_to Matrix{Integer} M
    return T.(transpose(Array(mat)))
end

## From Nemo matrices to Juila arrays
function convert_to_array(M::S, ::T=1) where {S <: MatElem,T}
    return T.(Array(M))
end

