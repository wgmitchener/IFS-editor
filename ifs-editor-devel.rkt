#! /usr/bin/racket
#lang racket/gui

; IFS Editor, November 2012 version
; Copyright 2012 by William Garrett Mitchener

;    This program is free software: you can redistribute it and/or modify
;    it under the terms of the GNU General Public License as published by
;    the Free Software Foundation, either version 3 of the License, or
;    (at your option) any later version.
;
;    This program is distributed in the hope that it will be useful,
;    but WITHOUT ANY WARRANTY; without even the implied warranty of
;    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
;    GNU General Public License for more details.
;
;    You should have received a copy of the GNU General Public License
;    along with this program.  If not, see <http://www.gnu.org/licenses/>.

(require racket/flonum)
(require racket/serialize)
(require framework)

(define vc-id "$Id: ifs-editor-devel.rkt,v dfbaa1d8cd01 2012-11-08 15:41:05Z garrett $")

(struct affine-transformation
  (x0 y0 a b c d)
  #:prefab)

(struct plane-point
  (x y)
  #:prefab)

(define (make-plane-point x0 y0)
  (plane-point
   (real->double-flonum x0)
   (real->double-flonum y0)))

(define (plane-point->complex p)
  (make-rectangular
   (real->double-flonum (plane-point-x p))
   (real->double-flonum (plane-point-y p))))

(define (complex->plane-point z)
  (make-plane-point (real-part z) (imag-part z)))

(define (plane-point-elementwise op)
  (lambda (p q)
    (plane-point
     (op (plane-point-x p) (plane-point-x q))
     (op (plane-point-y p) (plane-point-y q))) ))

(define plane-point-add (plane-point-elementwise fl+))
(define plane-point-subtract (plane-point-elementwise fl-))

(define (plane-point-scale alpha0 p)
  (define alpha (real->double-flonum alpha0))
  (plane-point
   (fl* alpha (plane-point-x p))
   (fl* alpha (plane-point-y p))))

(define (plane-point-dot u v)
  (fl+
   (fl* (plane-point-x u) (plane-point-x v))
   (fl* (plane-point-y u) (plane-point-y v))))

(define (plane-point-norm-sq u)
  (plane-point-dot u u))

(define (plane-point-cross u v)
  (fl-
   (fl* (plane-point-x u) (plane-point-y v))
   (fl* (plane-point-y u) (plane-point-x v))))
  
(define (apply-affine-transformation at pt)
  (define x0 (affine-transformation-x0 at))
  (define y0 (affine-transformation-y0 at))
  (define a (affine-transformation-a at))
  (define b (affine-transformation-b at))
  (define c (affine-transformation-c at))
  (define d (affine-transformation-d at))
  (define x (plane-point-x pt))
  (define y (plane-point-y pt))
  (plane-point 
   (fl+ x0 (fl* a x) (fl* b y))
   (fl+ y0 (fl* c x) (fl* d y))) )

(define (affine-transformation-fixed-point at)
  (define x0 (affine-transformation-x0 at))
  (define y0 (affine-transformation-y0 at))
  (define a (affine-transformation-a at))
  (define b (affine-transformation-b at))
  (define c (affine-transformation-c at))
  (define d (affine-transformation-d at))
  (define a1 (fl- a 1.0))
  (define d1 (fl- d 1.0))
  (define det1 (fl- (fl* a1 d1) (fl* b c)))
  (cond
    [(fl< det1 1e-8)
     #f]
    [else
     (plane-point
      (fl/ (fl- (fl* b y0) (fl* d1 x0)) det1)
      (fl/ (fl- (fl* c x0) (fl* a1 y0)) det1))]))

(define affine-transformation-identity 
  (affine-transformation 0.0 0.0 1.0 0.0 0.0 1.0))

(define (affine-transformation-compose at1 at2)
  (define x1 (affine-transformation-x0 at1))
  (define y1 (affine-transformation-y0 at1))
  (define a1 (affine-transformation-a at1))
  (define b1 (affine-transformation-b at1))
  (define c1 (affine-transformation-c at1))
  (define d1 (affine-transformation-d at1))

  (define x2 (affine-transformation-x0 at2))
  (define y2 (affine-transformation-y0 at2))
  (define a2 (affine-transformation-a at2))
  (define b2 (affine-transformation-b at2))
  (define c2 (affine-transformation-c at2))
  (define d2 (affine-transformation-d at2))

  (affine-transformation
   (fl+ x1 (fl+ (fl* a1 x2) (fl* b1 y2))) ;x3
   (fl+ y1 (fl+ (fl* c1 x2) (fl* d1 y2))) ;y3
   (fl+ (fl* a1 a2) (fl* b1 c2)) ;a3
   (fl+ (fl* a1 b2) (fl* b1 d2)) ;b3
   (fl+ (fl* c1 a2) (fl* d1 c2)) ;c3
   (fl+ (fl* c1 b2) (fl* d1 d2)) ;d3
   ) )

(define (determinant at)
  (define x0 (affine-transformation-x0 at))
  (define y0 (affine-transformation-y0 at))
  (define a (affine-transformation-a at))
  (define b (affine-transformation-b at))
  (define c (affine-transformation-c at))
  (define d (affine-transformation-d at))
  (fl- (fl* a d) (fl* b c)) )

(define (trace at)
  #;(define x0 (affine-transformation-x0 at))
  #;(define y0 (affine-transformation-y0 at))
  (define a (affine-transformation-a at))
  #;(define b (affine-transformation-b at))
  #;(define c (affine-transformation-c at))
  (define d (affine-transformation-d at))
  (fl+ a d))

(define (eigenvalues at)
  (define t (trace at))
  (define d (determinant at))
  (define disc (fl- (fl* t t) (fl* 4 d)))
  (list
   (fl/ (fl+ (fl- t) (sqrt disc)) 2)
   (fl/ (fl- (fl- t) (sqrt disc)) 2)) )

#;(define (singular-values at)
   (match at
    ((affine-transformation x0 y0 a b c d)
     (define A (+ (* a a) (* b b) (* c c) (* d d)))
     (define D (- (* a d) (* b c)))
     (define Q (sqrt (- (* A A) (* 4 D))))
     (list
      (sqrt (/ (+ A Q) 2))
      (sqrt (/ (- A Q) 2))))))

(define (orientation at)
  (let ([d (determinant at)])
    (if (fl< 0.0 d) -1
        (if (fl> 0.0 d) 1 0))))

(define (find-first predicate xs not-found)
  (if (null? xs) not-found
      (if (predicate (car xs))
          (car xs)
          (find-first predicate (cdr xs) not-found))))

(define drag-handler<%>
  (interface ()
    begin-drag
    continue-drag
    end-drag))

(define location<%>
  (interface ()
    get-location
    set-location))

(define size<%>
  (interface ()
    get-size
    set-size))

(define region<%>
  (interface ()
    contains-point?))

(define draw<%>
  (interface ()
    draw))

(define active<%>
  (interface ()
    activate
    deactivate))

(define observable<%>
  (interface ()
    add-observer
    remove-observer
    notify-observers))

; An observer is a function of type
; (object * hint-key * hint -> Any)
; hint-key is a symbol indicating the type of the event
(define observable-mixin
  (mixin () (observable<%>)
    
    (define observers (make-hasheq))
    
    (define/public (add-observer hint-keys observer)
      (for-each
       (lambda (hint-key)
         (hash-update! 
          observers hint-key
          (lambda (other-observers-box) 
            (box (cons observer (unbox other-observers-box))))
          ; If no observers yet, start with:
          (box '()) ))
       hint-keys) )
    
    (define/public (remove-observer observer)
      (hash-for-each
       observers
       (lambda (hint-key observers-box)
         (set-box!
          observers-box
          (remove observer (unbox observers-box))))) )
      
    (define/public (notify-observers hint-key hint)
      (define relevant-observers (unbox (hash-ref observers hint-key (box '()))))
      #;(display (format "notify-observers: ~s ~s ~s\n" this hint-key hint))
      #;(displayln relevant-observers)
      (for-each
       (lambda (observer)
         (observer this hint-key hint))
       relevant-observers))
    
    (super-new)))

(define drag-dispatcher%
  (class*
      (observable-mixin
       object%)
      
    (draw<%> active<%> region<%>)
    
    (init-field
     (managed-objects '())
     (active-object 'none))
    
    (inherit notify-observers)
    
    (define drag-in-progress-flag #f)

    (define/public (no-object-is-active?)
      (eq? 'none active-object))
    (define/public (an-object-is-active?)
      (not (no-object-is-active?)))
    
    (define (observable-callback object hint-key hint)
      (notify-observers 'subobject (list object hint-key hint)))
    
    (define (find-object-at x y)
      (find-first
        (lambda (obj) 
          (send obj contains-point? x y))
        managed-objects
        'none))
    
    (define/public (get-active-object)
      active-object)
    
    (define/public (contains-point? x y)
      (not (eq? 'none (find-object-at x y))))

    (define/public (clear)
      (deactivate 0 0)
      (define (loop)
        (cond
          [(not (eq? '() managed-objects))
           (remove-managed-object (car managed-objects))
           (loop)]))
      (loop))
    
    (define/public (add-managed-object object)
      (set! managed-objects
            (cons object managed-objects))
      #;(send object add-observer '(location) observable-callback)
      (notify-observers 'new-subobject object))
    
    (define/public (remove-managed-object object)
      (send object remove-observer observable-callback)
      (set! managed-objects
            (remq object managed-objects))
      (notify-observers 'remove-subobject object))
    
    (define/public (remove-active-object)
      (cond
        [(an-object-is-active?)
         (define o active-object)
         (define o-loc (send o get-location))
         (deactivate (plane-point-x o-loc) (plane-point-y o-loc))
         (remove-managed-object o)]))
    
    (define/public (dispatch-event mouse-event)
      (define is-dragging?
        (send mouse-event dragging?))
      (define x (send mouse-event get-x))
      (define y (send mouse-event get-y))
      (define button-change?
        (send mouse-event button-changed? 'left))
      (define click-begin? 
        (and button-change?
             (send mouse-event button-down? 'left)))
      (define click-end?
        (and button-change?
             (send mouse-event button-up? 'left)))
     
          
      (cond
        ; If the button is just clicked and no object is active, look for one
        [(and click-begin?
              (no-object-is-active?))
         (activate x y)]
        
        ; If the button is just clicked not on the active object, deactivate it
        ; This is a bit weird -- doesn't work if the active object is under
        ; another object and the click on the other object is also inside
        ; the active object
        [(and click-end?
              (an-object-is-active?)
              (not (send active-object contains-point? x y))
              (not (drag-in-progress?)))
         (deactivate x y)
         (activate x y)])
      
      (cond       
        ; If the button is pressed on the active object, begin a drag
        [(and click-begin?
              (an-object-is-active?)
              (send active-object contains-point? x y))
         (begin-drag x y)]
      
        ; If the button is down and we're dragging, continue drag:
        [(and is-dragging?
              (an-object-is-active?)
              (drag-in-progress?))
         (continue-drag x y)]
        
        ; If the button just came up and we're dragging, end drag:
        [(and click-end?
              (an-object-is-active?)
              (drag-in-progress?))
         (end-drag x y)]
        ))

    (define/public (activate x y)
      (define object-under-click (find-object-at x y))
 
      (cond
        [(not (eq? active-object object-under-click))
         (deactivate x y)
         (set! active-object object-under-click)
         (set! managed-objects
               (cons active-object
                     (remq active-object managed-objects)))])

      (cond
        [(an-object-is-active?)
          (send active-object activate x y)
          (send this notify-observers 'activate active-object)]))
    
    (define/public (deactivate x y)
      (cond
        [(an-object-is-active?)
          (send active-object deactivate x y)
          (send this notify-observers 'deactivate active-object)])
      (set! active-object 'none))
    
    (define/public (begin-drag x y)
      (cond
        [(an-object-is-active?)
         (set! drag-in-progress-flag #t)
         (send active-object begin-drag x y)
         (send this notify-observers 'begin-drag active-object)]))
    
    (define (drag-in-progress?)
      drag-in-progress-flag)
    
    (define/public (continue-drag x y)
      (cond
        [(an-object-is-active?)        
         (send active-object continue-drag x y)
         (send this notify-observers 'continue-drag active-object)]))
    
    (define/public (end-drag x y)
      (cond
        [(an-object-is-active?)
         (send active-object end-drag x y)
         (set! drag-in-progress-flag #f)
         (send this notify-observers 'end-drag active-object)]))
    
    (define/public (draw dc)
      (map
       (lambda (item)
         (send item draw dc))
       (reverse managed-objects))) ; ??? so I have a proper stack of managed objects
    
    (super-new) ))

(define has-location-mixin
  (mixin (observable<%>) (location<%>)
    (init-field 
     (location (make-plane-point 0 0)))
    (super-new)
    
    (inherit
      notify-observers)
    
    (define/public (set-location p)
      (define previous-location location)
      (set! location p)
      (notify-observers 'location previous-location))
    
    (define/public (get-location)
      location)
    ))

(define move-when-dragged-mixin
  (mixin (location<%>) (drag-handler<%>)
    
    (define previous-drag-point 'none)
    
    (define/public (begin-drag x y)
      (set! previous-drag-point 
            (make-plane-point x y)))
    
    (define/public (continue-drag x y)
      (define current-drag-point (make-plane-point x y))
      (define delta (plane-point-subtract current-drag-point previous-drag-point))
      (send this set-location
            (plane-point-add (send this get-location) delta))
      (set! previous-drag-point current-drag-point))
    
    (define/public (end-drag x y)
      (continue-drag x y)
      (set! previous-drag-point 'none)) 
    
    (super-new)))

(define handle-base%
  (class*
      (move-when-dragged-mixin 
       (has-location-mixin
        (observable-mixin object%)))
    
    (active<%>)
    
    (init-field
     (color "blue")
     (size 12)
     (active-state #f))
    
    (define/public (activate x y)
      (set! active-state #t))
    
    (define/public (deactivate x y)
      (set! active-state #f))
    
    (define/public (active?)
      active-state)
        
    (super-new) ))

(define square-handle%
  (class*
      handle-base%
    (draw<%> region<%>)
    
    (inherit-field
     color
     location
     size
     active-state)

    (define/public (contains-point? x y)
      (define x0 (plane-point-x location))
      (define y0 (plane-point-y location))
      (and
       (<= (- x0 size) x (+ x0 size))
       (<= (- y0 size) y (+ y0 size))) )

    (define/public (draw dc)
      (define x0 (plane-point-x location))
      (define y0 (plane-point-y location))
      (define p (new dc-path%))
      (send p rectangle (- x0 size) (- y0 size) (* 2 size) (* 2 size))
      
      (send* dc
        (set-brush color 'solid)
        (set-pen color 2 
                 (if active-state 'solid 'transparent))
        (draw-path p)))
    
    (super-new) ))

(define circle-handle%
  (class*
      handle-base%
    (draw<%> region<%>)
    
    (inherit-field
     location
     color
     size
     active-state)
    
     (define/public (contains-point? x y)
      (<= (plane-point-norm-sq
           (plane-point-subtract
            (make-plane-point x y)
            location))
          (* size size)))
      
    (define/public (draw dc)
      (define x0 (plane-point-x location))
      (define y0 (plane-point-y location))
      (define p (new dc-path%))
      (send p ellipse (- x0 size) (- y0 size) (* 2 size) (* 2 size))
      (send* dc
        (set-brush color 'solid)
        (set-pen color 2 
                 (if active-state 'solid 'transparent))
        (draw-path p)))
    
    (super-new) ))

(define circle-indicator%
  (class circle-handle%
    (define/override (contains-point? x y)
      #f)))

(define editable-state<%>
  (interface (active<%> drag-handler<%> draw<%> region<%>)))

(define inactive-editable-state%
  (class* object% (editable-state<%>)
    (init-field parent)
    (super-new)
    
    (define/public (activate x y)
      #f)
    (define/public (deactivate x y)
      #f)
    (define/public (begin-drag x y)
      #f)
    (define/public (continue-drag x y)
      #f)
    (define/public (end-drag x y)
      #f)
    (define/public (draw dc)
      #f)
    (define/public (contains-point? x y)
      #f) ))

(define active-editable-state%
  (class* drag-dispatcher% (editable-state<%>)
    (init-field
     parent
     (u-color "purple")
     (v-color "orange")
     #;(fp-color (make-object color% 0 80 0))
     )
    
    (inherit
      add-managed-object
      activate
      deactivate)
    
    (inherit-field
      active-object)

    (super-new)

    (define location-handle
      (new square-handle%
           [location (send parent get-corner-1)]))

    (define u-handle
      (new circle-handle% 
           [color u-color]
           [location (send parent get-corner-2)]))
    
    (define v-handle
      (new circle-handle%
           [color v-color]
           [location (send parent get-corner-4)]))
     
    (define (location-handle-observer handle location-hint-key previous-location)
      (send parent set-location (send handle get-location))
      ; that updates the corners, so:
      (send u-handle set-location (send parent get-corner-2))
      (send v-handle set-location (send parent get-corner-4)) )
    
    
    (define (u-handle-observer handle location-hint-key previous-location)
      (define u-loc (send u-handle get-location))
      (send parent set-corner-2 u-loc)
      (cond
        [(send u-handle active?)                    
         ; Also move v handle:
         (define origin (plane-point->complex
                         (send location-handle get-location)))
         (define old-z (- (plane-point->complex previous-location)
                          origin))
         (define new-z (- (plane-point->complex u-loc)
                          origin))
         (define old-w (- (plane-point->complex
                           (send v-handle get-location))
                          origin))
         (cond
           [(not (zero? new-z))
            (define new-w (* old-w (/ new-z old-z)))             
            #;(displayln (list origin old-z new-z old-w new-w))
            (send v-handle set-location 
                  (complex->plane-point (+ origin new-w))) ] )
         ]))
    
    (define (v-handle-observer handle hint-key previous-location)
      (send parent set-corner-4 (send v-handle get-location)))
    
    (send location-handle add-observer '(location) location-handle-observer)
    (send u-handle add-observer '(location) u-handle-observer)
    (send v-handle add-observer '(location) v-handle-observer)


    (define/override (begin-drag x y)
      (activate x y)
      (super begin-drag x y)
      (send parent notify-observers 'begin-sequence 'drag))
    
    (define/override (end-drag x y)
      (super end-drag x y)
      (deactivate x y)
      (send parent notify-observers 'end-sequence 'drag))
    
    (add-managed-object location-handle)
    (add-managed-object u-handle)
    (add-managed-object v-handle)
    ))

(define block%
  (class*
      (has-location-mixin
       (observable-mixin object%))
    
    (region<%> draw<%> active<%> drag-handler<%>)
    
    (define transparent-red (make-object color% 255 0 0 0.2))
    (define transparent-blue (make-object color% 0 0 255 0.2))
    
    (inherit-field location)
    (inherit
      notify-observers)
    
    (init-field
     (u (make-plane-point 100 0))
     (v (make-plane-point 0 -100))
     (u-color "purple")
     (v-color "orange"))
    
    (super-new)    
    
    (define inactive-state 
      (new inactive-editable-state%
           [parent this]))
    
    (define active-state
      (new active-editable-state%
           [parent this]
           [u-color u-color]
           [v-color v-color]))
    
    (define state inactive-state)
    
    (define/public (get-corner-1)
      location)
    
    (define/public (get-corner-2)
      (plane-point-add location u))
    
    (define/public (get-corner-3)
      (plane-point-add
       location 
       (plane-point-add u v)))
    
    (define/public (get-corner-4)
      (plane-point-add location v))
    
    (define/public (set-u p)
      (set! u p)
      (notify-observers 'u #f))
    
    (define/public (set-corner-2 p)
      (set-u (plane-point-subtract p location)))
    
    (define/public (set-v p)
      (set! v p)
      (notify-observers 'v p))
    
    (define/public (set-corner-4 p)
      (set-v (plane-point-subtract p location)))
    
    (define/public (basic-contains-point? x y)
      (define x0 (plane-point-x location))
      (define y0 (plane-point-y location))
      (define w (make-plane-point (- x x0) (- y y0)))
      (define u-component
        (/ (plane-point-dot u w)
           (plane-point-norm-sq u)))
      (define v-component
        (/ (plane-point-dot v w)
           (plane-point-norm-sq v)))
      (and
       (<= 0 u-component 1)
       (<= 0 v-component 1)))

      
    (define/public (contains-point? x y)
      (or (basic-contains-point? x y)
          (send state contains-point? x y)))
    
    (define/public (activate x y)
      (set! state active-state)
      (send state activate x y))
    
    (define/public (deactivate x y)
      (send state deactivate x y)
      (set! state inactive-state))

    (define/public (begin-drag x y)
      (send state begin-drag x y))
    
    (define/public (continue-drag x y)
      (send state continue-drag x y))
    
    (define/public (end-drag x y)
      (send state end-drag x y))
   
    (define/public (draw dc)
      (basic-draw dc)
      (send state draw dc))
    
    (define/public (basic-draw dc)
      (define location (send this get-location))
      (define a0 (plane-point-x location))
      (define a1 (plane-point-y location))
      (define u0 (plane-point-x u))
      (define u1 (plane-point-y u))
      (define v0 (plane-point-x v))
      (define v1 (plane-point-y v))
 
      (define main-outline (new dc-path%))
      (send* main-outline
        (move-to a0 a1)
        (line-to (+ a0 u0) (+ a1 u1))
        (line-to (+ a0 u0 v0) (+ a1 u1 v1))
        (line-to (+ a0 v0) (+ a1 v1))
        (close))
      
      (define u-edge (new dc-path%))
      (send* u-edge
        (move-to a0 a1)
        (line-to (+ a0 u0) (+ a1 u1)))

      (define v-edge (new dc-path%))
      (send* v-edge
        (move-to a0 a1)
        (line-to (+ a0 v0) (+ a1 v1)))
      
      ; Fill color is red for positive orientation,
      ; blue for negative orientation.
      ; Screen coordinates have a negative orientation,
      ; so this is flipped:
      (define color
        (if (> (plane-point-cross u v) 0)
            transparent-blue
            transparent-red))
      
      (send* dc
        (set-brush color 'solid)
        (set-pen color 2 'solid)
        (draw-path main-outline)
        (set-brush color 'transparent)
        (set-pen u-color 2 'solid)
        (draw-path u-edge)
        (set-pen v-color 2 'solid)
        (draw-path v-edge) ))
    
 ))


(define managed-canvas%
  (class canvas%
    
    (init-field
     (unit-size 300)
     (padding 150)
     (u-color "purple")
     (v-color "orange")
     (dispatcher (new drag-dispatcher%)))
    
    (inherit
      refresh)
    
    (define/override (on-event mouse-event)
      (send dispatcher dispatch-event mouse-event))
    
    (super-new
     (min-width (+ padding unit-size padding))
     (min-height (+ padding unit-size padding))
     (paint-callback
          (lambda (this dc)
            (define main-outline (new dc-path%))
            (define u+padding (+ padding unit-size))
            ; The unit square
            (send* main-outline
              (move-to padding padding)
              (line-to padding u+padding)
              (line-to u+padding u+padding)
              (line-to u+padding padding)
              (close))
            
            (define u-edge (new dc-path%))
            (send* u-edge
              (move-to padding u+padding)
              (line-to u+padding u+padding))

            (define v-edge (new dc-path%))
            (send* v-edge
              (move-to padding u+padding)
              (line-to padding padding))
      

            (send* dc 
              (set-smoothing 'smoothed)
              (set-pen "black" 3 'solid)
              (set-brush "black" 'transparent)
              (draw-path main-outline)
              (set-pen u-color 3 'solid)
              (draw-path u-edge)
              (set-pen v-color 3 'solid)
              (draw-path v-edge))
            
            (send dispatcher draw dc))))
    
    (define (do-refresh object hint-key hint)
      (refresh))
    
    (send dispatcher add-observer 
          '(subobject
            new-subobject
            remove-subobject
            activate
            deactivate
            begin-drag
            continue-drag
            end-drag) do-refresh)
    ))

; west: screen horizontal coordinate of the left side of the unit box
; south: screen vertical coordinate of the bottom side of the unit box
; unit-size: screen pixel size of the unit box
(define (block->affine-transformation west-0 south-0 unit-size-0 block)
  (define n real->double-flonum)
  (define west (n west-0))
  (define south (n south-0))
  (define unit-size (n unit-size-0))
  ; screen-pixel coordinates:
  (define U (get-field u block))
  (define V (get-field v block))
  (define Q (send block get-location))
  (define Qx (n (plane-point-x Q)))
  (define Qy (n (plane-point-y Q)))

  ; Cartesian plane coordinates:
  (define q
    (plane-point
     (fl/ (fl- Qx west) unit-size)
     (fl/ (fl- Qy south) unit-size)))
  (define u
    (plane-point-scale (fl/ 1.0 unit-size) U))
  (define v
    (plane-point-scale (fl/ 1.0 unit-size) V))
  (affine-transformation
   (n (plane-point-x q))
   (n (- (plane-point-y q)))
   (n (plane-point-x u))
   (n (plane-point-x v))
   (n (- (plane-point-y u)))
   (n (- (plane-point-y v)))) )

(define (blocks->IFS west south unit-size blocks)
  (map
   (lambda (b)
     (block->affine-transformation west south unit-size b))
   blocks))

(define ifs-canvas%
  (class canvas%
    (init-field
     picture-size
     size-threshold
     (color "black")
     (max-depth 20)
     (affine-transformations '())
     (when-done-callback (lambda () #f))
     )
    
    (inherit
      refresh)
    
    (define bitmap #f)
    
    (define/public (run-when-done-callback)
      (when-done-callback))
     
    (define/public (set-affine-transformations at)
      (set! affine-transformations at)
      (set! bitmap #f)
      (fork-render))

    (define (fork-render)
      (thread
       (lambda ()
         (set! bitmap 
           (render-fractal-to-bitmap
            affine-transformations
            picture-size
            (/ picture-size 2)
            max-depth
            (real->double-flonum size-threshold)))
         (refresh)
         (run-when-done-callback)
         )))
 
    (super-new
     [min-width (* 2 picture-size)]
     [min-height (* 2 picture-size)]
     [paint-callback
        (lambda (c dc)
          (cond
            [bitmap
             (send dc draw-bitmap bitmap 0 0) ]))])

    (fork-render)
    
    ))

(define (draw-ifs dc at ats depth size-threshold)
  (cond
    ; If it's small, just draw the box:
    [(or (= 0 depth)
         (fl< (abs (determinant at)) size-threshold))
     (draw-box-for dc at)]
    [else
     (for-each
      (lambda (t)
        (draw-ifs 
         dc 
         (affine-transformation-compose at t)
         ats (- depth 1) size-threshold))
      ats)] ))

(define (render-fractal-to-bitmap 
         affine-transformations
         picture-size
         padding
         max-depth
         size-threshold)

  (define color "black")
  (define size (inexact->exact 
                (ceiling (+ (* 2 padding) picture-size))))
  (define bitmap 
    (make-object bitmap% size size #t))
  
  (define dc
    (new bitmap-dc% [bitmap bitmap]))
  
  (send* dc
    (translate padding (+ padding picture-size))
    (scale picture-size (- picture-size))
    (set-pen color 0 'transparent)
    (set-brush color 'solid))
  
  (draw-ifs dc affine-transformation-identity affine-transformations
            max-depth (real->double-flonum size-threshold))
  
  bitmap)

(define (draw-box-for dc at)
  (define p (new dc-path%))
  (match at
    [(affine-transformation x0 y0 a b c d)
     (send* p
       (move-to x0 y0)
       (line-to (+ x0 a) (+ y0 c))
       (line-to (+ x0 a b) (+ y0 c d))
       (line-to (+ x0 b) (+ y0 d))
       (close))
     (send dc draw-path p)]))

(define-syntax-rule
  (number-editor parent-object field-label variable)
  
  (new text-field%
       [parent parent-object]
       [label field-label]
       [min-width 200]
       [init-value (format "~s" variable)]
       [callback
        (lambda (tf event)
          (define v (send tf get-value))
          (define n (string->number v))
          (cond
            [n
             (set! variable (real->double-flonum n))]
            ))]) )

(define (run-render-dialog parent ats)
  
  (define picture-size 256)
  (define padding 128)
  (define max-depth 24)
  
  (define df
    (new dialog%
         [parent parent]
         [min-width 300]
         [label "Render fractal"]))
  
  (define parameter-v-box
    (new vertical-panel%
         [parent df]))
  
  
  (number-editor parameter-v-box "Unit square pixel size" picture-size)
  (number-editor parameter-v-box "Padding pixels" padding)
  (number-editor parameter-v-box "Maximum recursion depth" max-depth)        
  
  (define bottom-box
    (new horizontal-panel%
         [parent parameter-v-box]))
  
  (new button%
       [parent bottom-box]
       [label "Cancel"]
       [callback
        (lambda (button event)
          (send df show #f))])
  
  (new button%
       [label "Render"]
       [parent bottom-box]
       [callback
        (lambda (button event)
          (define file-path 
            (finder:put-file "Untitled.png" (current-directory)))
          (send df set-cursor (make-object cursor% 'watch))
          (define bitmap 
            (render-fractal-to-bitmap ats picture-size padding max-depth
                                      (/ (* picture-size picture-size))))
          (send bitmap save-file file-path 'png)
          (send df show #f)
          )])
  
  (send df show #t)
  )

(define (run-affine-transformation-dialog parent at)
  (define x0 (affine-transformation-x0 at))
  (define y0 (affine-transformation-y0 at))
  (define a (affine-transformation-a at))
  (define b (affine-transformation-b at))
  (define c (affine-transformation-c at))
  (define d (affine-transformation-d at))

  (define result #f)
  
  (define df
    (new dialog%
         [parent parent]
         [label "Edit affine transformation"]))
  
  (define v-box
    (new vertical-panel%
         [parent df]))

  (number-editor v-box "x0" x0)
  (number-editor v-box "y0" y0)
  (number-editor v-box "a" a)
  (number-editor v-box "b" b)
  (number-editor v-box "c" c)
  (number-editor v-box "d" d)
  
  (define bottom-box
    (new horizontal-panel%
         [parent v-box]))
  
  (new button%
       [parent bottom-box]
       [label "Cancel"]
       [callback
        (lambda (button event)
          (send df show #f))])
  
  (new button%
       [parent bottom-box]
       [label "Set"]
       [callback
        (lambda (button event)
          (set! result
                (affine-transformation x0 y0 a b c d))
          (send df show #f))])
  
  (send df show #t)
  
  result)

(define (run-ifs-editor)
  (define unit-size 256)
  (define padding 128)
  (define west padding)
  (define half-unit (/ unit-size 2))
  (define quarter-unit (/ unit-size 4))
  (define south (+ padding unit-size))
  (define picture-size unit-size)
  (define drawing-size-threshold (/ (* unit-size unit-size)))
  (define picture-max-depth 24)
  (define preview-picture-size (/ unit-size 2))
  (define preview-size-threshold (* 32 drawing-size-threshold))
  (define ats '())
    
  (define (new-default-block m)
    (send* m 
      (add-managed-object
       (new block% 
            [location (make-plane-point 
                       (+ west quarter-unit)
                       (- south quarter-unit))]
            [u (make-plane-point half-unit 0)]
            [v (make-plane-point 0 (- half-unit))]))
      (activate (+ west half-unit) (- south half-unit))))
 
  (define (affine-transformation->block at)
    (define x0 (affine-transformation-x0 at))
    (define y0 (affine-transformation-y0 at))
    (define a (affine-transformation-a at))
    (define b (affine-transformation-b at))
    (define c (affine-transformation-c at))
    (define d (affine-transformation-d at))
    
    (new block%
         [location
          (plane-point-add
           (make-plane-point west south)
           (plane-point-scale unit-size (plane-point x0 (- y0))))]
         [u (plane-point-scale unit-size (plane-point a (- c)))]
         [v (plane-point-scale unit-size (plane-point b (- d)))]))
  
  (define (sierpinski-blocks m)
    (send* m
      (add-managed-object
       (new block%
            [location (make-plane-point west south)]
            [u (make-plane-point half-unit 0)]
            [v (make-plane-point 0 (- half-unit))]))
      (add-managed-object
       (new block%
            [location (make-plane-point
                       (+ west half-unit)
                       south)]
            [u (make-plane-point half-unit 0)]
            [v (make-plane-point 0 (- half-unit))]))
      (add-managed-object
       (new block%
            [location (make-plane-point west (- south half-unit))]
            [u (make-plane-point half-unit 0)]
            [v (make-plane-point 0 (- half-unit))])) ))
      
  
  (define frame
    (new frame%
         [label "Canvas"]))

  (define v-box
    (new vertical-panel%
         [parent frame]))
  
  (define button-box
    (new horizontal-panel%
         [parent v-box]))

  (define h-box-main
    (new horizontal-panel%
         [parent v-box]))
  
  (define c
    (new managed-canvas%
         [unit-size unit-size]
         [padding padding]
         [stretchable-height #f]
         [stretchable-width #f]
         [parent h-box-main]))
  
  (define side-bar-v-box
    (new vertical-panel%
         [parent h-box-main]))
  
  (define t
    (new text-field%
         [label #f]
         [parent side-bar-v-box]
         [enabled #t]
         [min-width 300]
         [min-height 150]
         [style '(multiple)]))
  
  (define edit-button
    (new button%
         [label "Edit"]
         [parent side-bar-v-box]
         [enabled #f]
         [callback
          (lambda (button event)
            (define block (send m get-active-object))
            (cond
              [(not (eq? 'none block))
               (define new-at 
                 (run-affine-transformation-dialog
                  frame
                  (block->affine-transformation west south unit-size block)))
               (cond
                 [new-at
                  (send m remove-active-object)
                  (send m add-managed-object
                        (affine-transformation->block new-at))]) ]))] ))
  
  (define preview-canvas
    (new ifs-canvas%
         [parent side-bar-v-box]
         [picture-size preview-picture-size]
         [size-threshold preview-size-threshold]
         [color (make-object color% 0 0 0 0.2)]
         [max-depth 3]
         [affine-transformations ats]))
  
  (define m (get-field dispatcher c))
  
  (define (build-ifs)
    (blocks->IFS west south unit-size (get-field managed-objects m)))

  (define (update-ifs-detail-display block hint-key hint)
    (define object (send m get-active-object))
    (define at (block->affine-transformation west south unit-size object))
    (define fp (affine-transformation-fixed-point at))
    (define x0 (affine-transformation-x0 at))
    (define y0 (affine-transformation-y0 at))
    (define a (affine-transformation-a at))
    (define b (affine-transformation-b at))
    (define c (affine-transformation-c at))
    (define d (affine-transformation-d at))

    (send t set-value
          (string-append
           (format 
            "x0: ~a\ny0: ~a\na: ~a\nb: ~a\nc: ~a\nd: ~a\ndet: ~a\n"
            x0 y0 a b c d
            (determinant at))
           (if fp
               (format "fp: (~a,~a)\n"
                       (plane-point-x fp)
                       (plane-point-y fp) )
               ""))) )
          
  (define (update-preview object hint-key hint)
    (send preview-canvas set-affine-transformations (build-ifs)))
  
  (send m add-observer '(new-subobject) update-preview)
  (send m add-observer '(remove-subobject) update-preview)
  
  (send m add-observer '(activate)
        (lambda (m activate-hint-key block)
          (send edit-button enable #t)
          (send block add-observer '(end-sequence) update-ifs-detail-display)
          (send block add-observer '(end-sequence) update-preview)
          (update-ifs-detail-display block 'u #f) ))
  
  (send m add-observer '(deactivate)
        (lambda (m deactivate-hint-key block)
          (send edit-button enable #f)
          (send t set-value "")
          (send block remove-observer update-ifs-detail-display)
          (send block remove-observer update-preview)
          #f))
  
  ; to get started:
  (sierpinski-blocks m)

  (new button%
       [parent button-box]
       [label "Clear"]
       [callback
        (lambda (button event)
          (send m clear))])
  
  (new button%
       [parent button-box]
       [label "Load"]
       [callback
        (lambda (button event)
          (define file-name 
            (finder:get-file (current-directory)))
          (cond
            [file-name
             (define ats-mementos
               (with-input-from-file file-name read))
             (send m clear)
             (for-each
              (lambda (at-m)
                (define at (deserialize at-m))
                (send m add-managed-object
                      (affine-transformation->block at)))
              ats-mementos)]) )])

    (new button%
       [parent button-box]
       [label "Save"]
       [callback
        (lambda (button event)
          (define file-name 
            (finder:put-file "Untitled.ifs" (current-directory)))
          (with-output-to-file file-name
            (lambda ()
              (write (map serialize (build-ifs)))
              )
            #:exists 'replace)
          )])

  (new button%
       [parent button-box]
       [label "Add"]
       [callback
        (lambda (button event)
          (new-default-block m))])
  
  (new button%
       [parent button-box]
       [label "Delete"]
       [callback
        (lambda (button event)
          (send m remove-active-object))])
  
  (new button%
       [parent button-box]
       [label "Sketch"]
       [callback
        (lambda (button event)
          (define f
            (new frame% [label "Fractal"]))
          (send f set-cursor (make-object cursor% 'watch))
          
          (define c
            (new ifs-canvas%
                 [parent f]
                 [picture-size picture-size]
                 [size-threshold drawing-size-threshold]
                 [max-depth picture-max-depth]
                 [affine-transformations (build-ifs)]
                 [when-done-callback
                  (lambda () (send f set-cursor #f))])) 
          (send f show #t))])

  (new button%
       [parent button-box]
       [label "Render"]
       [callback
        (lambda (button event)
          (run-render-dialog frame (build-ifs)))])
           
  (send frame show #t))

; main

(run-ifs-editor)
