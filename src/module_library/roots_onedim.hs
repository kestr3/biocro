-- householder methods
newton :: (Fractional t) => (t -> t) -> (t -> t) -> t -> t
newton fun jac x = x - fun x / jac x

halley :: (Fractional a) => (a -> a) -> (a -> a) -> (a -> a) -> a -> a
halley f df df2 x = x - r * recip (1 - r * s / 2)
  where
    r = f x / df x
    s = df2 x / df x

steffensen :: (Fractional a) => (a -> a) -> a -> a
steffensen fun x = x - y / dy
  where
    y = fun x
    dy = fun (x + y) / y - 1

-- multi step
twoStep :: ((a -> a) -> a -> a -> a) -> (a -> a) -> (a, a) -> (a, a)
twoStep formula fun (x0, x1) = (x1, formula fun x0 x1)

secantFormula :: (Fractional a) => (a -> a) -> a -> a -> a
secantFormula fun x0 x1 = x1 - y1 / slope
  where
    slope = (y1 - y0) / (x1 - x0)
    y1 = fun x1
    y0 = fun x0

secant :: (Fractional a) => (a -> a) -> (a, a) -> (a, a)
secant fun (x0, x1) = twoStep secantFormula fun (x0, x1)

threeStep :: ((a -> a) -> a -> a -> a -> a) -> (a -> a) -> (a, a, a) -> (a, a, a)
threeStep formula fun (x0, x1, x2) = (x1, x2, formula fun x0 x1 x2)

-- lagrange polynomial cenetered at x=a
lagrangeQuad :: (Fractional a) => a -> a -> a -> a -> a
lagrangeQuad a b c x = (x - b) * (x - c) / (a - b) / (a - c)

-- derived using Lagrange interpolation
invQuadInterpFormula :: (Fractional a) => (a -> a) -> a -> a -> a -> a
invQuadInterpFormula fun x0 x1 x2 = w0 * x0 + w1 * x1 + w2 * x2
  where
    w0 = lagrangeQuad y0 y1 y2 0
    w1 = lagrangeQuad y1 y0 y2 0
    w2 = lagrangeQuad y2 y0 y1 0
    y0 = fun x0
    y1 = fun x1
    y2 = fun x2

inverseQuadraticInterpolation :: (Fractional a) => (a -> a) -> (a, a, a) -> (a, a, a)
inverseQuadraticInterpolation = threeStep invQuadInterpFormula

-- bracketing methods
bisection :: (Fractional b, Ord a, Num a) => (b -> a) -> (b, b) -> (b, b)
bisection fun (left, right) = updateBracket fun c (left, right)
  where
    c = (left + right) / 2

regulaFalsi :: (Fractional b, Ord b) => (b -> b) -> (b, b) -> (b, b)
regulaFalsi fun (left, right) = updateBracket fun c (left, right)
  where
    c = secantFormula fun left right

dekker fun (b', a, b) = (b, a1, b1)
  where
    (a1, b1) = if abs (fun ap) < abs (fun bp) then (bp, ap) else (ap, bp)
    ap = if not $ sameSign (fun b) (fun bp) then b else a
    bp = if liesBetween s (m, b) then s else m
    m = (a + b) / 2
    s = secantFormula fun b' b

liesBetween :: (Ord a) => a -> (a, a) -> Bool
liesBetween s (a, b) =
  if a < b then (a < s) && (s < b) else (b < s) && (s < a)

updateBracket fun x (left, right) =
  if sameSign (fun x) (fun left)
    then (x, right)
    else (left, x)

sameSign :: (Ord a, Num a) => a -> a -> Bool
sameSign x y = x * y > 0

isZero :: (Ord a, Num a) => a -> a -> Bool
isZero tol x = abs x < tol

isClose :: (Ord a, Num a) => a -> a -> a -> a -> Bool
isClose atol rtol x y = abs (x - y) < atol + rtol * max x y

-- root at x = 1
testFun :: (Fractional a, Ord a) => a -> a
testFun x = if x > 0 then (x ^ 4 - 1) / 4 else -(1 / 4)

testFunDeriv :: (Num a, Ord a) => a -> a
testFunDeriv x = if x > 0 then x ^ 3 else 0

testFun2ndDeriv :: (Num a, Ord a) => a -> a
testFun2ndDeriv x = if x > 0 then 3 * x ^ 2 else 0
