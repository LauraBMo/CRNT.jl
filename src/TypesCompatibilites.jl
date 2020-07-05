
## Add float arithmetic compatibility with AbstractAlgebra rational type
Base.Float64(x::Nemo.fmpq) = Base.Float64(Rational(x))
Base.Float32(x::Nemo.fmpq) = Base.Float32(Rational(x))
Base.Float16(x::Nemo.fmpq) = Base.Float16(Rational(x))
Base.BigFloat(x::Nemo.fmpq) = Base.BigFloat(Rational(x))

## In case we want to push! a fmpq coeff to a float array
Base.convert(::Type{T}, x::Nemo.fmpq) where {T<:AbstractFloat} = T(Rational(x))

## From float to fmpq
(a::Nemo.FlintRationalField)(x::T) where {T<:AbstractFloat} = Nemo.fmpq(rationalize(x))
Nemo.fmpq(x::T) where {T<:AbstractFloat} = Nemo.fmpq(rationalize(x))
