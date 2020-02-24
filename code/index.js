
import { Renderer } from "./renderer.js";
import { get, on } from "./aliases.js";
import { Ship } from "./entity/ship.js";
import { Planetoid } from "./entity/planetoid.js";
import { ShipController } from "./controller/shipcontroller.js";
import { InputManager } from "./input.js";
import API from "./api.js";

let api = new API();
api.p2 = p2;
api.pworld = new p2.World({
  gravity: [0, 0],
  applyDamping:false,
  applyGravity:false
});
api.canvas = get("canvas");
api.input = new InputManager(document.body);
api.renderer = new Renderer(api, true, 1);

for (let i=0; i<200; i++) {
let planet1 = new Planetoid(api);
planet1.x = Math.random()*500;
planet1.y = Math.random()*500;
api.renderer.addEntity(planet1);
}

let ship = new Ship(api);
console.log(ship);
let shipController = new ShipController(api, ship);
api.renderer.addEntity(ship);

on(api.canvas, "wheel", (evt) => {
  evt.preventDefault();
  api.renderer.addZoom((evt.deltaY * api.renderer.zoom) / 50);
});

let handleInput = () => {
  shipController.step();
  api.renderer.setCenter(ship.x, ship.y);
  api.renderer.setCursor(api.input.mouse.x, api.input.mouse.y);
};

setInterval(() => {
  handleInput();
}, 1000 / 60);
