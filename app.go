package main

import "fmt"
import "os"
import "bufio"
import "regexp"
import "math"

type crystal struct { a, b, c, α, β, γ float64 }

var cs = map[int]string {
    1: "Triclinic",
    2: "Monoclinic",
    3: "Orthorhombic",
    4: "Tetragonal",
    5: "Trigonal",
    6: "Hexagonal",
    7: "Cubic",
}
var prompt string = "> "
var pi float64 = math.Pi

func rad(in float64) (float64) { return (in * pi/180.0) }
func deg(in float64) (float64) { return (in * 180.0/pi) }

// takes deg. as in/out, computes in rads:

func sin(in float64) (float64)  { return deg(math.Sin(rad(in))) }
func cos(in float64) (float64)  { return deg(math.Cos(rad(in))) }
func tan(in float64) (float64)  { return deg(math.Tan(rad(in))) }
func isin(in float64) (float64) { return deg(math.Asin(rad(in))) }
func icos(in float64) (float64) { return deg(math.Acos(rad(in))) }
func itan(in float64) (float64) { return deg(math.Atan(rad(in))) }
func pow(b,e float64) (float64) { return math.Pow(b,e) }

// get stdin input, error check

func getLatticeParam(q string) (f float64) {
	stdin := bufio.NewReader(os.Stdin)
    fmt.Print(q+prompt)
    for {
    	_, err := fmt.Fscan(stdin, &f)
		if err == nil {	return f }
		stdin.ReadString('\n')
		fmt.Print(q+prompt)
	}
}

// get stdin input, error check + regular expressions

func getCrystalSys() (i int) {

    var s string
    var m bool
    stdin := bufio.NewReader(os.Stdin)

    fmt.Print(prompt)

    for {
        _, err := fmt.Fscan(stdin, &s)
        if err != nil {
            stdin.ReadString('\n')
            fmt.Print(prompt)
        }

        for i, w := range cs {
            if len(s) > 3 { m, _ = regexp.MatchString("(?i)"+w[0:4], s[0:4])
                } else { m = false }
            if m || 48 + i == int(s[0]) { return i }
        }

        stdin.ReadString('\n')
        fmt.Print(prompt)
    }
}

// program starts here!

func main() {

    var c crystal
    var h, k, l, ds, th float64

    fmt.Println("------ Please selet your crystal/lattice system ----")
    for i := 1; i < 8; i++ {
        fmt.Printf("   %d. "+cs[i]+"\n",i)
    }
    CrystalSys := getCrystalSys()
    fmt.Println("")

    fmt.Println("------ Please enter wavelength of x-rays -----------")
    Wavelength := getLatticeParam("λ (Å)\n")
    fmt.Println("")

    fmt.Println("------ Please enter your lattice parameters --------")
    switch CrystalSys {
        case 1: //tric
            c.a = getLatticeParam("a (Å) ≠ b, c\n")
            c.b = getLatticeParam("b (Å) ≠ a, c\n")
            c.c = getLatticeParam("c (Å) ≠ a, b\n")
            c.α = getLatticeParam("α (≠90 deg.)\n")
            c.β = getLatticeParam("β (≠90 deg.)\n")
            c.γ = getLatticeParam("γ (≠90 deg.)\n")
        case 2: //mono
            c.a = getLatticeParam("a (Å) ≠ b, c\n")
            c.b = getLatticeParam("b (Å) ≠ a, c\n")
            c.c = getLatticeParam("c (Å) ≠ a, b\n")
            c.α = 90.0
            c.β = getLatticeParam("β (≠90 deg.)\n")
            c.γ = 90.0
        case 3: //orth
            c.a = getLatticeParam("a (Å) ≠ b, c\n")
            c.b = getLatticeParam("b (Å) ≠ a, c\n")
            c.c = getLatticeParam("c (Å) ≠ a, b\n")
            c.α = 90.0
            c.β = 90.0
            c.γ = 90.0
        case 4: //tetr
            c.a = getLatticeParam("a, b (Å)    \n")
            c.b = c.a
            c.c = getLatticeParam("c (Å) ≠ a, b\n")
            c.α = 90.0
            c.β = 90.0
            c.γ = 90.0
        case 5: //trig
            c.a = getLatticeParam("a, b, c (Å) \n")
            c.b = c.a
            c.c = c.b
            c.α = getLatticeParam("α (≠90 deg.)\n")
            c.β = getLatticeParam("β (≠90 deg.)\n")
            c.γ = getLatticeParam("γ (≠90 deg.)\n")
        case 6: //hexa
            c.a = getLatticeParam("a, b (Å)    \n")
            c.b = c.a
            c.c = getLatticeParam("c (Å) ≠ a, b\n")
            c.α = 90.0
            c.β = 90.0
            c.γ = 120.0
        case 7: //cubi
            c.a = getLatticeParam("a, b, c (Å) \n")
            c.b = c.a
            c.c = c.b
            c.α = 90.0
            c.β = 90.0
            c.γ = 90.0
        default: return
    }
    fmt.Println("")

    fmt.Println("------ variables for computation -------------------")
    fmt.Printf("   %s\n",cs[CrystalSys])
    fmt.Printf("   a: %g Å\n",c.a)
    fmt.Printf("   b: %g Å\n",c.b)
    fmt.Printf("   c: %g Å\n",c.c)
    fmt.Printf("   α: %g deg.\n",c.α)
    fmt.Printf("   β: %g deg.\n",c.β)
    fmt.Printf("   γ: %g deg.\n",c.γ)
    fmt.Printf("   λ: %g Å\n",Wavelength)
    fmt.Println("")
    fmt.Println("------ h, k, l, d-spacing, 2theta ------------------")

    // d-spacing formulas from Dr. Falak Sher
    // National Workshop on Crystal Structure Determination using Powder XRD,
    // organized by the Khwarzimic Science Society, 15 – 17 August 2007 (http://www.khwarzimic.org).

    for l = 0.0; l < 3.0; l=l+1.0 {
        for k = 0.0; k < 3.0; k=k+1.0 {
            for h = 0.0; h < 3.0; h=h+1.0 {
                switch CrystalSys {
                    case 1: //tri
                        var vl float64 = math.Sqrt(1 - (pow(cos(c.α),2)) - (pow(cos(c.β),2)) - (pow(cos(c.γ),2)) + (2.0*cos(c.α)*cos(c.β)*cos(c.γ)))
                        var t1 float64 = 1/pow((c.a * c.b * c.c * vl),2)
                        var t2 float64 = pow(c.b,2) * pow(c.c,2) * pow(sin(c.α),2)
                        var t3 float64 = pow(c.a,2) * pow(c.c,2) * pow(sin(c.β),2)
                        var t4 float64 = pow(c.a,2) * pow(c.b,2) * pow(sin(c.γ),2)
                        var t5 float64 = pow(c.c,2) * c.a * c.b * (cos(c.α) * cos(c.β) - cos(c.γ))
                        var t6 float64 = pow(c.a,2) * c.c * c.b * (cos(c.β) * cos(c.γ) - cos(c.α))
                        var t7 float64 = pow(c.b,2) * c.a * c.c * (cos(c.γ) * cos(c.α) - cos(c.β))
                        ds = t1*(t2 * pow(h,2) + t3 * pow(k,2) + t4 * pow(l,2) + 2*t5*h*k + 2*t6*k*l + 2*t7*h*l)
                    case 2: //MONOCLINIC!
                        var t1 float64 = 1.0/pow(sin(c.β),2)
                        var t2 float64 = pow(h,2)/pow(c.a,2)
                        var t3 float64 = pow(k,2) * pow(sin(c.β),2) / pow(c.b,2)
                        var t4 float64 = pow(l,2)/pow(c.c,2)
                        var t5 float64 = (2.0 * h * l * cos(c.β)) / (c.a * c.c)
                        ds = t1 * (t2 + t3 + t4 - t5)
                    case 3: //ORTHORHOMBIC!
                        ds = (pow(h,2)/pow(c.a,2)) + (pow(k,2)/pow(c.b,2)) + (pow(l,2)/pow(c.c,2))
                    case 4: //TETRAGONAL!
                        ds = (pow(h,2) + pow(k,2))/pow(c.a,2) + pow(l,2)/pow(c.c,2)
                    case 5: //TRIGONAL!
                        var t1 float64 = (pow(h,2) + pow(k,2) + pow(l,2)) * pow(sin(c.α),2)
                        var t2 float64 = 2.0 * (h*k + k*l + h*l) * pow(cos(c.α),2)
                        var t3 float64 = cos(c.α)
                        var t4 float64 = pow(c.a,2) * (1.0 - 3.0 * pow(cos(c.α),2) + 2.0 * pow(cos(c.α),3))
                        ds = (t1 + t2 - t3) / t4
                    case 6: //HEXAGONAL!
                        ds = 4.0/3.0 * ((pow(h,2) + h*k + pow(k,2))/pow(c.a,2)) + pow(l,2)/pow(c.c,2)
                    case 7: //CUBIC!
                        ds = ( pow(h,2) + pow(k,2) + pow(l,2) ) / pow(c.a,2)
                    default: return
                }
                ds = 1.0/math.Sqrt(ds)
                th = deg(isin(Wavelength/(2.0 * ds)))
                fmt.Printf("%.1g %.1g %.1g   %1.4f    %1.4f\n",h,k,l,ds,th*2.0)
            }
        }
    }
}