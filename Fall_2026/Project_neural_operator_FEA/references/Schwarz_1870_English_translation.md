# On a Limiting Process by an Alternating Method

**H. A. Schwarz**

*(From a lecture given on 30 May.)*

Vierteljahrsschrift der Naturforschenden Gesellschaft in Zürich, vol. 15 (1870), pp. 272–286.
Original title: "Ueber einen Grenzübergang durch alternirendes Verfahren."

English translation from the Internet Archive / Biodiversity Heritage Library scan (item `vierteljahrsschr15natu`, PDF pages 280–294). The original is in the public domain. Translator's notes are in square brackets. Where the scan's text recognition garbled a formula, the formula has been reconstructed from the mathematics and is marked [reconstructed].

---

The mode of inference known by the name of Dirichlet's principle, which in a certain sense must be regarded as the foundation of the branch of the theory of analytic functions developed by Riemann, is subject, as is now generally admitted, to very well-founded objections with respect to rigour, and their complete removal has, as far as I know, not yet been achieved by the efforts of mathematicians.

By continuing some investigations concerning certain kinds of mapping problems, part of which has been published in volume 70 of Borchardt's Journal and in the essay "Zur Theorie der Abbildung" accompanying the programme of the Federal Polytechnic School for the winter semester 1869–70, I have been led to a method of proof by which, as I believe I have convinced myself, all theorems whose proof Riemann sought to carry out in his published papers by means of Dirichlet's principle can be proved rigorously.

The following communication is essentially an extract from a paper concerning the integration of the partial differential equation $\Delta u = 0$, which I communicated in November of last year to Mr. Kronecker and to some other mathematicians close to me.

The matter is essentially only to prove the existence of a function $u$ which, for a certain given region $T$ of the independent real variables $x$ and $y$, satisfies the partial differential equation

$$\Delta u = \frac{\partial^2 u}{\partial x^2} + \frac{\partial^2 u}{\partial y^2} = 0$$

and in addition certain prescribed boundary and discontinuity conditions.

For brevity I restrict myself here to the case in which the side conditions are boundary conditions only; that is, in which it is required that the function $u$ always be finite and have, along the boundary of the region $T$, prescribed finite values belonging to one or more continuous sequences. For the general case can be reduced to this case by the method to be communicated in what follows.

For the applicability of the method of proof in question it is by no means necessary to make the assumption that the boundary line of $T$ has only a finite number of corners and, in general, a finite definite radius of curvature at every point, an assumption which Messrs. Weber and Neumann have made in their investigations directed at the same goal. (See Borchardt's Journal, vol. 71, p. 29, and the Reports of the Mathematical-Physical Class of the Royal Saxon Society of Sciences, session of 21 April 1870.) Not even continuity in the change of direction of the tangent to the boundary line is required; it suffices rather to know that the boundary line can be divided into a finite number of pieces such that, in the interior of each of these pieces, the change of direction of the tangent always takes place in the same sense, even if it does so infinitely often by jumps, so that the boundary line may possess infinitely many corners.

Cusps of the boundary line are not excluded either. For those cusps which arise from the contact of two analytic lines having, in the neighbourhood of the point of contact, the character of algebraic curves, I have carried out the investigation; however, in order to avoid unnecessary lengthiness here, no account is taken of the occurrence of cusps in what follows.

The success of the proof, whose basic idea is communicated here, rests in the last instance on the following lemma:

> Let the boundary line of a region $T$, for which it is possible to integrate the partial differential equation $\Delta u = 0$ in accordance with arbitrary boundary conditions, be divided into a finite number of segments (parts). Let these be arranged in two groups, so that each group contains at least one segment. To the individual segments assign odd or even ordinal numbers, according as they belong to the first or the second group, and denote by $P$ the points which separate the segments with even ordinal number from those with odd ordinal number. In the interior of $T$ imagine a finite number of analytic lines $L$ given, which have in common with the segments of odd ordinal number either no point or only end points $P$ of these segments, without however touching them at these points.
>
> Now imagine a function $u$ determined for the region $T$ which satisfies the partial differential equation $\Delta u = 0$ and has, at every point of the boundary of $T$, the value $0$ or $+1$ according as the ordinal number of the segment in whose interior the point in question lies is even or odd. Then the upper limit, or respectively the maximum, of all the values which the function $u$ takes along the lines $L$ is a positive number $q$ which is smaller than $1$.
>
> If now, for the same region $T$, with the same division of the boundary into segments of even and odd ordinal number and for the same lines $L$, a function $u_1$ is determined which satisfies the differential equation $\Delta u_1 = 0$, has the value zero on the boundary of $T$ along the segments of even ordinal number, and whose value, arbitrarily prescribed along the segments of odd ordinal number, does not exceed the quantity $g$ in absolute value, then the absolute value of the values which the function $u_1$ can take at the points of the lines $L$ nowhere exceeds the value $g \cdot q$, where $q$ has the meaning given above and is therefore smaller than $1$.

For the area of a circle, and for all simply connected surfaces whose conformal mapping onto the area of a circle is known, the integration of the partial differential equation $\Delta u = 0$ under prescribed boundary conditions involves no difficulty. With regard to this problem, permission may be granted to refer to a small essay communicated in this journal (pp. 113–128 of the current volume); in it, interruptions of continuity in the sequence of values prescribed for the function $u$ along the boundary were, for brevity, expressly excluded; nevertheless the conclusions developed there hold, mutatis mutandis, even when the sequence of prescribed boundary values suffers an interruption of continuity at a finite number of boundary points.

Having shown that for a number of simpler regions the differential equation $\Delta u = 0$ can be integrated in accordance with arbitrary boundary conditions, the task is to prove that, also for a less simple region composed of those in a certain way, the integration of the differential equation in accordance with arbitrary boundary conditions is possible. For the proof of this theorem a limiting process can serve which has a great analogy with the well-known procedure for producing a rarefied space by means of a two-cylinder air pump. For in the one case as in the other, the period of the operation consists of two individual operations which come into effect alternately, and which have indeed the same purpose but are not identical with respect to the manner of their action, being rather symmetric in a certain sense.

Such a limiting process may briefly be called a **limiting process by an alternating method**.

Let two regions $T_1$ and $T_2$ be given which have one or more regions $T^*$ in common and whose boundary lines do not touch one another. (In the accompanying schematic figure $T_1$ is the area of a circle, $T_2$ the area of a square.) [The figure shows a disc and a square overlapping.]

Let the system of all parts of the boundary of $T_1$ which lie outside $T_2$ be denoted by $L_0$, and the system of all the remaining parts, lying inside $T_2$, by $L_2$.

Likewise the boundary of $T_2$ falls into the systems $L_1$ and $L_3$, namely if $L_1$ denotes the system of all pieces which lie inside the region $T_1$, and $L_3$ the system of all pieces which lie outside $T_1$.

It is assumed that both for the region $T_1$ and for the region $T_2$ it is possible to integrate the differential equation $\Delta u = 0$ in accordance with arbitrary boundary conditions; the task is to show that this is also possible for the region $T_1 + T_2 - T^* = T$, which contains the regions $T_1$ and $T_2$ as parts, but in which the region $T^*$ common to the regions $T_1$ and $T_2$ is to be counted only once.

Both for the region $T_1$ and the system $L_1$, and for the region $T_2$ and the system $L_2$, the conditions of the lemma mentioned above are fulfilled; in the first case let the system $L_2$, in the second the system $L_1$, take the place of the group of segments with even ordinal number. It is therefore possible to determine two numbers $q_1$ and $q_2$ which play the role of the number $q$ in the lemma and which are both smaller than $1$.

Keeping the above analogy, the region $T^*$ corresponds to the receiver of the air pump, the regions $T_1 - T^*$ and $T_2 - T^*$ to the interiors of the two pump cylinders, and the lines $L_1$ and $L_2$ to the valves.

Let the values of the function $u$ be prescribed arbitrarily on the boundary of $T$, that is along $L_0$ and $L_3$; let $g$ be the upper and $k$ the lower limit of these values; let the difference $g - k$ be denoted by $G$.

Now assume along $L_2$ an arbitrary sequence of values, for example the value $k$ at all points of $L_2$, and determine for the region $T_1$ a function $u_1$ which has the prescribed values along $L_0$, the value $k$ along $L_2$, and satisfies the differential equation $\Delta u_1 = 0$ in the interior of $T_1$. By the assumption made about the region $T_1$ there is such a function. (First stroke of the first piston.)

Imagine the values which the function $u_1$ has along $L_1$ fixed, and determine for the region $T_2$ a function $u_2$ which has the prescribed values along $L_3$, agrees along $L_1$ with the previously determined function $u_1$, and for which $\Delta u_2 = 0$. By the assumption made about the region $T_2$ there is such a function. (First stroke of the second piston.)

The value of $u_2 - u_1$, or of $u_2 - k$, along $L_2$ is smaller than $g - k = G$.

Now determine for the region $T_1$ a function $u_3$ which has the prescribed values along $L_0$, agrees along $L_2$ with $u_2$, and for which $\Delta u_3 = 0$. (Second stroke of the first piston.)

The difference $u_3 - u_1$ is at no point in the interior of $T_1$ negative, and is in absolute value smaller than $G$; along $L_1$, however, by the lemma mentioned, it is smaller than $G \cdot q_1$, because $u_3 - u_1$ has the value zero along $L_0$ and is smaller than $G$ along $L_2$.

Imagine the value of the function $u_3$ along $L_1$ fixed, and determine for the region $T_2$ a function $u_4$ which agrees along $L_1$ with $u_3$, has the prescribed values along $L_3$, and for which $\Delta u_4 = 0$. (Second stroke of the second piston.)

The difference $u_4 - u_2$ has the value zero along $L_3$, and along $L_1$, where it agrees with $u_3 - u_1$, it is positive and smaller than $G \cdot q_1$; therefore in the interior of $T_2$, $u_4 - u_2$ is nowhere negative and always smaller than $G \cdot q_1$, but along $L_2$ it is smaller than $G \cdot q_1 \cdot q_2$.

By continuing this alternating method one arrives at a sequence of infinitely many functions with odd and with even index. The former are defined for the region $T_1$, the latter for the region $T_2$, in such a way that they have the prescribed values along $L_0$ and $L_3$ respectively, and satisfy the partial differential equation $\Delta u = 0$ in the interior of the regions for which they are defined.

For the region $T^*$ both the functions with odd index and those with even index are defined, and indeed they agree with one another alternately along $L_1$ and along $L_2$. Namely, along $L_1$, $u_{2n-1} = u_{2n}$, and along $L_2$, $u_{2n+1} = u_{2n}$.

It is now not difficult to prove that the functions with odd index and those with even index approach without limit, as the index increases, definite limit functions $u'$ and $u''$, which are defined by the equations

$$u' = u_1 + (u_3 - u_1) + (u_5 - u_3) + \cdots + (u_{2n+1} - u_{2n-1}) + \cdots \text{ ad inf.}$$

$$u'' = u_2 + (u_4 - u_2) + (u_6 - u_4) + \cdots + (u_{2n+2} - u_{2n}) + \cdots \text{ ad inf.}$$

The series on the right-hand sides converge absolutely and, for all pairs of values $x, y$ under consideration, uniformly; for

$$(u_{2n+1} - u_{2n-1}) < G \cdot (q_1 q_2)^{\,n-1} \quad\text{and}\quad (u_{2n+2} - u_{2n}) < G \cdot (q_1 q_2)^{\,n-1} \cdot q_1 .$$

Both along $L_1$ and along $L_2$, $u' = u''$. In the interior of $T_1$, $\Delta u' = 0$; in the interior of $T_2$, $\Delta u'' = 0$; therefore for every point of $T^*$, $u' = u''$, because on the whole boundary of $T^*$ the two functions agree with one another.

The two functions $u'$ and $u''$ are therefore values of the same function $u$, which is defined for the whole region $T = T_1 + T_2 - T^*$, satisfies the partial differential equation $\Delta u = 0$ in its interior, and takes the prescribed values on the boundary $L_0 + L_3$.

With this the proof of the assertion stated above has been indicated: under the stated assumptions it is also possible, for the region $T$, to integrate the partial differential equation $\Delta u = 0$ in accordance with arbitrarily prescribed boundary conditions.

By repeated application and suitable modification of the limiting process by an alternating method just mentioned, the existence of a function $u$ for a given region can be established even when, besides the boundary conditions, discontinuity conditions are also prescribed, or, as in the case of the Abelian integrals, discontinuity conditions alone, in those cases for which Riemann asserted existence in his papers and sought to prove it by means of Dirichlet's principle.

The method of proof explained extends not only to the case in which the simply or multiply connected Riemann surface representing the region $T$ is contained, in its whole extent, in the same plane or on the same spherical surface, but holds essentially unchanged also when this surface is spread out over a polyhedral surface composed of several plane or spherical faces.

With the help of this extension one can, among other things, also carry out the proof that a simply connected region spread out over such a polyhedral surface can be mapped conformally onto the area of a circle if it possesses a boundary line returning into itself, and onto the surface of a sphere if the region is a simply connected and closed region.

With this the question of the possibility of the determination of the constants, to which the conformal mapping of the simply connected surface of a polyhedron bounded by plane faces onto a sphere can be reduced (see Borchardt's Journal, vol. 70, p. 119), is also answered.

A special case of the mapping problem just mentioned occurs when it is a matter of mapping conformally the simply connected area of a plane polygon bounded by straight segments onto the area of a circle, whether the area of the polygon lies entirely in the finite part of the plane or contains the point at infinity of the plane once or several times in its interior; branch points in its interior are not excluded either. In this problem the only difficulty likewise lies in proving the possibility of determining a certain number of constants, partly real and partly conjugate complex, on which the mapping function is made to depend, in such a way that all the conditions of the problem are satisfied.

This difficulty can be overcome by applying a method developed by Mr. Weierstrass. The application of the above limiting method offers a new means of overcoming it.

The same holds for the proof of the possibility of that determination of constants to which the problem of the conformal mapping of a simply connected figure bounded by arcs of circles onto the area of a circle is reduced (loc. cit., p. 117). In this case too the surface may contain branch points or the point at infinity in its interior.

**Postscript.** Recently three papers of Mr. Christoffel have become known to me (Annali di Matematica, edited by Brioschi and Cremona, vol. IV, pp. 1–9; Nachrichten von der K. Gesellschaft der Wissenschaften zu Göttingen 1870, pp. 283–298 and 359–369), whose content, in so far as it relates to the above communication, gives occasion for some remarks.

On p. 1 of the fourth volume of the Annali one reads: "… the determination of the stationary temperatures on a rectangular surface $F$ offers no difficulty, the associated problem of the stationary temperatures on the complementary surface $F'$ has remained completely inaccessible to the methods employed so far …"

and on p. 284 of the Nachrichten mentioned: "… one then arrives at a remarkable class of problems which presents such striking difficulties that the solution of a problem of this kind has hitherto, in spite of all efforts, succeeded only in the single, entirely elementary case in which the boundary of $\mathfrak{F}_1$ (the surface extending to infinity in all directions which is left of a plane when a simply connected finite piece $\mathfrak{F}$ is cut out of it) is a circle. In particular, as will find its explanation in what follows, all attempts have failed to treat any one of the special, especially interesting cases in which $\mathfrak{F}_1$ is bounded by a rectilinear figure."

In the face of this assertion it seems appropriate to draw attention to some simple examples, which perhaps also offer some interest in themselves to those to whom they are new.

1. In the plane of the complex quantity $z$ let a parabola be given whose focus is the point $z = 0$ and whose vertex is the point $z = +1$. By the function $Z = \tan^2\!\left(\tfrac{\pi}{4}\sqrt{z}\right)$ the interior, and by the function $Z = \dfrac{2}{\sqrt{z}} - 1$ the exterior, of the area of the parabola is mapped connectedly and similarly in the smallest parts [that is, conformally] onto the area of a circle, by which, as is well known, the corresponding heat problem for the exterior and for the interior of the parabola is to be regarded as solved. [Both formulas reconstructed; the scan is illegible at these two places. With $w = \sqrt{z}$ the parabola $y^2 = 4(1-x)$ is the image of the lines $\operatorname{Re} w = \pm 1$, so the interior corresponds to the strip $|\operatorname{Re} w| < 1$, which $\tan(\pi w/4)$ carries to the unit disc, and the exterior to the half-plane $\operatorname{Re} w > 1$, which $2/w - 1$ carries to the unit disc.]

2. Let the equation of an ellipse $\dfrac{x^2}{a^2} + \dfrac{y^2}{b^2} = 1$ be given, with $a^2 - b^2 = 1$. By the function $Z = \operatorname{sin\,am}\!\left(\dfrac{2K}{\pi}\arcsin z\right)$, with nome $q = \left(\dfrac{a-b}{a+b}\right)^2$, the interior, and by the function $Z = \dfrac{z - \sqrt{z^2 - 1}}{a - b}$ the exterior, of the ellipse is mapped conformally onto the area of a circle, by which the corresponding heat problem is to be regarded as solved for this case too. [The arguments of the elliptic function and of the exterior map are reconstructed; "sin am" is Jacobi's notation for the elliptic sine, sn.]

3. Let a square be given whose corners are the four points $z = +1, +i, -1, -i$. By the function $Z = \operatorname{sin\,am}(Kz)$, $k = \tfrac{1}{\sqrt 2}$, the interior of the square is mapped conformally onto the area of a circle [modulus reconstructed], and by the function

$$z = \frac{1}{C}\int \frac{\sqrt{1 + Z_1^4}}{Z_1^2}\, dZ_1, \qquad C = \int_0^{\pi/2} \sqrt{\sin\varphi}\; d\varphi$$

conversely the area of the circle of radius 1 described about the origin in the plane of the complex quantity $Z_1$ is mapped conformally onto the exterior of that square, if the complex variable of integration $Z_1$ is restricted to that circular area and the constant of integration is determined from the condition that $\lim\,(z - \tfrac{1}{C Z_1})$ for $Z_1 = 0$ is also equal to zero. [The integrand and the normalising condition are reconstructed; the scan shows only "$C = \int \sqrt{\sin\varphi}\,d\varphi$" and "lim $(z - \ldots)$ for $Z_1 = 0$" clearly.] By the mappings given, however, the heat problem for the interior and for the exterior of the square may be regarded as solved. (Compare Borchardt's Journal, vol. 70, p. 115.)

4. Let a circular-arc triangle be given. The conformal mapping of the area of a circle onto the interior and onto the exterior of a circular-arc triangle can be carried out without difficulty by means of hypergeometric series. (Loc. cit., p. 117.)

The number of these examples could be increased further.

It must not be left unmentioned here that considerable objections are also to be raised against the remaining content of the papers named.

Of a complete proof one would doubtless have to demand that it also be proved that it is in fact possible to determine all the constants on whose determination the matter depends in such a way that all the conditions of the problem are satisfied.

It is further to be noted that in the two papers mentioned in the Göttinger Nachrichten the investigation is expressly restricted to the case in which the boundary line of the region to be mapped is given by an "indecomposable equation". By this, therefore, both the cases in which the boundary line consists of several pieces of different analytic lines, and the case in which this line at no place has the character of an algebraic curve, are excluded from the outset.

Finally it is not to be overlooked that in a very large number of cases, to which the simple case of the conformal mapping of the interior of an ellipse onto the interior of a circle also belongs, the final formulas set up in the papers named remain burdened with an unresolved difficulty, a circumstance which makes it appear questionable whether these formulas may be regarded as a proof of the possibility of solving the general problem posed.

August 1870.

---

### Translator's notes

- **The lemma** is the whole engine. In modern language: a harmonic function with boundary data $1$ on the "far" boundary and $0$ on the interface has maximum $q < 1$ on any interior curve kept away from the "far" boundary; by linearity and the maximum principle, every subdomain solve therefore contracts the interface error by a fixed factor. The product $q_1 q_2$ is the contraction factor of one full alternating cycle, exactly the quantity that governs the convergence of every Schwarz method since.
- **The air pump.** Schwarz's analogy is a two-cylinder vacuum pump: the overlap $T^*$ is the receiver (bell jar), the non-overlapping parts of $T_1$ and $T_2$ are the two cylinders, and the interfaces $L_1$ and $L_2$ are the valves. Each "stroke" is one subdomain solve.
- **Setting.** Schwarz's goal was not a numerical method but a rigorous existence proof for the Dirichlet problem of Laplace's equation on domains assembled from pieces (a disc and a square in his figure), to replace Riemann's use of Dirichlet's principle, which Weierstrass had shown to be unjustified. The method became a computational tool only in the 1980s, after P.-L. Lions reinterpreted it.
- **Borchardt's Journal** is the Journal für die reine und angewandte Mathematik (Crelle's Journal), edited by Borchardt at the time. "Sin am" is Jacobi's amplitude notation for the elliptic function now written $\operatorname{sn}$; $K$ is the complete elliptic integral of the first kind.
- **Companion paper.** The essay on pp. 113–128 of the same volume, "Ueber die Integration der partiellen Differentialgleichung $\Delta u = 0$ für die Fläche eines Kreises", supplies the disc solution (the Poisson integral) that the alternating method builds on.
