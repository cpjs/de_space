package main

import "fmt"
import "os"
import "bufio"
import "regexp"
import "math"

type crystal struct {
	a, b, c, α, β, γ 	float64
	symmetry 			string
}

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

func main() {

    var c crystal
    var h, k, l, ds float64

    fmt.Println("------ Please selet your crystal/lattice system ----")
    for i := 1; i < 8; i++ {
        fmt.Printf("   %d. "+cs[i]+"\n",i)
    }
    CrystalSys := getCrystalSys()
    fmt.Println("")

    fmt.Println("------ Please enter wavelength of x-rays -----------")
    Wavelength := getLatticeParam("λ (nm)\n")
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
    fmt.Printf("   λ: %g nm\n",Wavelength)
    fmt.Println("")
    fmt.Println("------ h, k, l, d-spacing, 2theta ------------------")

    for l = 0.0; l < 3.0; l=l+1.0 {
        for k = 0.0; k < 3.0; k=k+1.0 {
            for h = 0.0; h < 3.0; h=h+1.0 {
                switch CrystalSys {
                    case 1: //tri
                        var v1 float64 = math.Pow(math.Cos(c.α),2)
                        var v2 float64 = math.Pow(math.Cos(c.β),2)
                        var v3 float64 = math.Pow(math.Cos(c.γ),2)
                        var f4 float64 = 2*math.Cos(c.α)*math.Cos(c.β)*math.Cos(c.γ)
                        var v5 float64 = math.Sqrt(1 - v1 - v2 - v3 + v4)
                        var t1 float64 = 1/math.Pow((c.a*c.b*c.c*v5),2)
                        var t2 float64 = math.Pow(c.b,2)*math.Pow(c.c,2)*math.Pow(math.Sin(c.α),2)
                        var t3 float64 = math.Pow(c.a,2)*math.Pow(c.c,2)*math.Pow(math.Sin(c.β),2)
                        var t4 float64 = math.Pow(c.a,2)*math.Pow(c.b,2)*math.Pow(math.Sin(c.γ),2)
                        var t5 float64 = math.Pow(c.c,2)*c.a*c.b*(math.Cos(c.α)*math.Cos(c.β)-math.Cos(c.γ))
                        var t6 float64 = math.Pow(c.a,2)*c.c*c.b*(math.Cos(c.β)*math.Cos(c.γ)-math.Cos(c.α))
                        var t7 float64 = math.Pow(c.b,2)*c.a*c.c*(math.Cos(c.γ)*math.Cos(c.α)-math.Cos(c.β))
                        ds = t1*(t2*math.Pow(h,2) + t3*math.Pow(k,2) + t4*math.Pow(l,2) + 2*t5*h*k + 2*t6*k*l + 2*t7*h*l)
                    case 2: // MONOCLINIC
                        var t1 float64 = 1/(math.Pow(math.Sin(c.β),2))
                        var t2 float64 = math.Pow(h,2)/math.Pow(c.a,2)
                        var t3 float64 = (math.Pow(k,2)*math.Pow(math.Sin(c.β),2))/(math.Pow(c.b,2))
                        var t4 float64 = math.Pow(l,2)/math.Pow(c.c,2)
                        var t5 float64 = (2*h*l*math.Cos(c.β))/(c.a*c.c)
                        ds = t1*(t2 + t3 + t4 - t5)
                    case 3: // ORTHORHOMBIC
                        var t1 float64 = math.Pow(h,2)/math.Pow(c.a,2)
                        var t2 float64 = math.Pow(k,2)/math.Pow(c.b,2)
                        var t3 float64 = math.Pow(l,2)/math.Pow(c.c,2)
                        ds = t1 + t2 + t3
                    case 4: //TETRAGONAL
                        var t1 float64 = (math.Pow(h,2)+math.Pow(k,2))/math.Pow(c.a,2)
                        var t2 float64 = math.Pow(l,2)/math.Pow(c.c,2)
                        ds = t1 + t2
                    case 5: //TRIGONAL
                        var t1 float64 = (math.Pow(h,2)+math.Pow(k,2)+math.Pow(l,2))*math.Pow(math.Sin(c.α),2)
                        var t2 float64 = 2.0*((h*k)+(k*l)+(h*l))*math.Pow(math.Cos(c.α),2)
                        var t3 float64 = math.Cos(c.α)
                        var t4 float64 = math.Pow(c.a,2)*(1-(3*math.Pow(math.Cos(c.α),2))+(2*math.Pow(math.Cos(c.α),3)))
                        ds = (t1 + t2 - t3)/t4
                    case 6: //HEXAGONAL
                        var t1 float64 = (4.0/3.0)*((math.Pow(h,2)+h*k+math.Pow(k,2))/math.Pow(c.a,2))
                        var t2 float64 = math.Pow(l,2)/math.Pow(c.c,2)
                        ds = t1 + t2
                    case 7: //cubic
                        ds = (math.Pow(h,2) + math.Pow(k,2) + math.Pow(l,2))/math.Pow(c.a,2)
                    default: return
                }
                ds = math.Sqrt(math.Pow(ds,-1))
                th2 := math.Pow(math.Sin(Wavelength/(0.2*ds)),-1)
                fmt.Printf("%.1g %.1g %.1g   %1.4f    %1.4f\n",h,k,l,ds,th2)
            }
        }
    }
}