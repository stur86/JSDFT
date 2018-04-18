// A style variable needs to be exposed, yet read only for safety

function get_style() {
    return $("#current_style").html();
}

function style_set(index) {

    var old_style = get_style();
    var new_style = $($(".style_select").find(".button")[index]).html();

    $("#current_style").html(new_style);
    $("body").removeClass("style_" + old_style);
    $("body").addClass("style_" + new_style);

    // Edit the various labels
    $(".particles_name").html(new_style);

    // Also re-init everything
    clear_all();

}