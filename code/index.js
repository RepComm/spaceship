
import { Renderer } from "./renderer.js";
import { get, on } from "./aliases.js";
import { Ship } from "./entity/ship.js";
import { Planetoid } from "./entity/planetoid.js";
import { ShipController } from "./controller/shipcontroller.js";
import { InputManager } from "./input.js";
import { Circle } from "./colliders.js";

let canvas = get("canvas");
let inputManager = new InputManager(document.body);
let renderer = new Renderer(canvas, true, 1);

let ship = new Ship();
ship.collider = new Circle();
ship.collider.radius = 0.2;
let shipController = new ShipController(inputManager, ship);

console.log(ship);

for (let i=0; i<100; i++) {
    let pl = new Planetoid();
    pl.collider = new Circle();
    pl.x = (Math.random()*200)*2-1;
    pl.y = (Math.random()*200)*2-1;
    pl.mass = Math.random()*10;
    pl.color = "#" + parseInt(Math.random()*0xffff);
    pl.scale = pl.mass/5;
    pl.collider.radius = pl.scale;
    pl.velocity.x = Math.random()*2-1;
    pl.velocity.y = Math.random()*2-1;
    pl.rotationalVelocity = Math.random()/8;
    renderer.addEntity(pl);
}

renderer.addEntity(ship);

on(canvas, "wheel", (evt) => {
    evt.preventDefault();
    renderer.addZoom((evt.deltaY * renderer.zoom) / 50);
});

let handleInput = () => {
    shipController.step();
    renderer.setCenter(ship.x, ship.y);
    renderer.setCursor(inputManager.mouse.x, inputManager.mouse.y);
};

setInterval(() => {
    handleInput();
}, 1000 / 60);
