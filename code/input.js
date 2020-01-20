
import { on } from "./aliases.js";

class InputManager {
    /**
     * @param {HTMLElement} element 
     */
    constructor (element) {
        this.element = element;

        /**@type {Map<String, Boolean>} */
        this.keys = new Map();
        on(this.element, "keyup", (e)=>this.onKeyUp(e));
        on(this.element, "keydown", (e)=>this.onKeyDown(e));

        this.mouse = {
            x:0,
            y:0
        };
        this.mouse.buttons = new Map();
        on(this.element, "mousemove", (e)=>this.onMouseMove(e));
        on(this.element, "mouseup", (e)=>this.onMouseUp(e));
        on(this.element, "mousedown", (e)=>this.onMouseDown(e));
    }

    onKeyUp (evt) {
        this.keys.set(evt.key, false);
    }

    onKeyDown (evt) {
        this.keys.set(evt.key, true);
    }

    onMouseMove (evt) {
        this.mouse.x = evt.layerX;
        this.mouse.y = evt.layerY;
    }

    /**
     * @param {MouseEvent} evt 
     */
    onMouseUp (evt) {
        this.mouse.buttons.set(evt.button, false);
    }
    /**
     * @param {MouseEvent} evt 
     */
    onMouseDown (evt) {
        this.mouse.buttons.set(evt.button, true);
    }

    /**Is a mouse button pressed down?
     * @param {Number} button
     * @returns {Boolean}
     */
    isMousePressed (button=0) {
        return this.mouse.buttons.get(button) || false;
    }

    /**Is a key pressed down?
     * @param {String} name
     * @returns {Boolean}
     */
    isKeyPressed (name) {
        return this.keys.get(name) || false;
    }
}

export { InputManager };