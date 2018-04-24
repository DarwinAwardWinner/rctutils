# Inverse of sitools::f2si
si2f <- function(string, unit="") {
    if (length(string) == 0) {
        return(numeric(0))
    }
    sifactor <- c(1e-24, 1e-21, 1e-18, 1e-15, 1e-12, 1e-09, 1e-06,
                  0.001, 1, 1000, 1e+06, 1e+09, 1e+12, 1e+15, 1e+18, 1e+21,
                  1e+24)
    pre <- c("y", "z", "a", "f", "p", "n", "u", "m",
             "", "k", "M", "G", "T", "P", "E", "Z", "Y")

    rx <- rex(
        ## Leading whitespace
        start,
        zero_or_more(space),

        ## Capture a floating point number
        capture(
            ## Sign
            maybe(one_of("+", "-")),
            ## Integer part
            zero_or_more(digit),
            ## Decimal point
            maybe("."),
            ## Fractional part (or integer part when decimal is not
            ## present)
            one_or_more(digit),
            ## Exponential notation
            maybe(
                one_of("e", "E"),
                maybe(one_of("+", "-")),
                one_or_more(digit)
            )
        ),

        ## Space between number and unit
        zero_or_more(space),

        ## Capture SI prefix
        capture(maybe(one_of(pre))),

        unit,

        ## Trailing whitespace
        zero_or_more(space),
        end
    )

    m <- str_match(string, rx)
    base <- as.numeric(m[,2])
    p <- m[,3]
    fac <- sifactor[match(p, pre)]
    base * fac
}
