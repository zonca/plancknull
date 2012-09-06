$(document).ready(function(){
    //remove right click menu
    $("a.thumbnail").bind("contextmenu", function(e) {
    return false;
    });

    $("#clearstorage").click(function() {
        $.each($.jStorage.index(), function(index, key) { 
                $.jStorage.deleteKey(key);
            });
            createSelected();
    });

    //$(".nav a").click(function() {
    //    $(".active").removeClass("active");
    //    $(this).parent("li.dropdown").addClass("active");
    //}) 

    function createSelected()
    {
        var output = $('#output');
        output.empty();
        var $selectedArray = $.jStorage.get("selected", new Array())
        $(".clicked").removeClass("clicked");
        $.each($.jStorage.index(), function(index, key) { 
            var $imgclone = $($.jStorage.get(key))
            output.append($imgclone);
            $(".table #" + key).addClass('clicked');
        });
        $("a[rel^='prettyPhoto[Sel]']").prettyPhoto({'social_tools':''});
        $("a[rel^='prettyPhoto[Sel]']").bind("contextmenu", function(e) {
        return false;
        });
        $("a[rel^='prettyPhoto[Sel]']").mousedown(function(e){
            if (e.which == 3) {
                $.jStorage.deleteKey($(this).attr('id').slice(0, -4));
                createSelected();
            }
        });
    }

    createSelected();

    $('.mapstable a.thumbnail').mousedown(function(e){
        if (e.which == 3) {
            var $img = $(this).toggleClass('clicked');
            var $imgId = $img.attr('id');
            var $value = $.jStorage.get($imgId, false)
            if($value)
            {
                $.jStorage.deleteKey($imgId);
            }
            else
            {
                var $imgclone = $img.clone();
                $imgclone = $imgclone.removeClass('thumbnail');
                $imgclone = $imgclone.attr('id', $imgId + '_sel');
                $imgclone = $imgclone.removeClass('clicked');
                $imgclone = $imgclone.attr('rel', 'prettyPhoto[Sel]');
                $.jStorage.set($imgId, $imgclone.clone().wrap('<p>').parent().html());
            }
            createSelected();
        }
    });

    function colorSummaryTable()
    {
        // color summary table
        $('.summarytable td').each(function() {
            $(this).css('background-color', 'white');
            var fValue = parseFloat($(this).text());
            var orangeValue = parseFloat($('#orange').val());
            if (fValue > orangeValue) {
                $(this).css('background-color', 'orange');
            };
            var redValue = parseFloat($('#red').val());
            if (fValue > redValue) {
                $(this).css('background-color', 'red');
            };
        });
    };
    colorSummaryTable();
    $('.updatesummarycolor').change(colorSummaryTable);
});
