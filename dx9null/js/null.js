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

    var output = $('#output');

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

    $('a.thumbnail').mousedown(function(e){
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
});
