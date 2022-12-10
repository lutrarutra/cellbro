
// function make_resizable() {

//     const resizer = document.querySelector('.resizer');
//     const top = document.querySelector('.top');
//     const bottom = document.querySelector('.bottom');
//     bottom.innerHTML = "";
//     const main = document.querySelector('.main');
//     const secondary = document.querySelector('.secondary');

//     const main_plot = document.querySelector('.main-plot');
//     const secondary_plot = document.querySelector('.secondary-plot');
//     const bottom_plot = document.querySelector('.bottom-plot');

//     const main_figure = document.querySelector('.main-figure');
//     const secondary_figure = document.querySelector('.secondary-figure');
//     const bottom_figure = document.querySelector('.bottom-figure');

//     const body = document.querySelector('body');

//     var w = window.innerWidth;
//     var h = window.innerHeight;
//     var x = 0;
//     var y = 0;
//     var clicked = false;

    

//     resizer.addEventListener("mousedown", (e) => {
//         clicked = true;
//         x = e.clientX;
//         y = e.clientY;
//         window.addEventListener("mousemove", mouse_move);
//     });

//     window.addEventListener("mouseup", (e) =>{
//         if (!clicked) return;
//         window.removeEventListener("mousemove", mouse_move);
//         var top_rect = top.getBoundingClientRect();

//         console.log(e.clientY);
//         console.log(top_rect.bottom);
//         console.log(top_rect.height);

//         top.style.height = e.clientY + "px";
//         clicked = false;
//     });

//     function mouse_move(e) {
//         resizer.style.position = "absolute";
//         resizer.style.left = e.clientX  + "px";
//         resizer.style.top = e.clientY  + "px";
//     }

// }


// window.addEventListener("load", (event) => {
//     setTimeout(make_resizable, 1000);
// });

